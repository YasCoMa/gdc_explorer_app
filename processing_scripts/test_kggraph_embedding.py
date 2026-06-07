import os
import re
import json
import uuid
import asyncio
import requests
import datetime
import numpy as np
import pandas as pd
import networkx as nx

from gql import Client, gql
from graphql import print_schema
from gql.transport.aiohttp import AIOHTTPTransport
from subprocess import run, PIPE
from tqdm import tqdm


import torch
import torch_geometric
device = 'cuda' if torch.cuda.is_available() else 'cpu'
torch.serialization.add_safe_globals([torch_geometric.data.hetero_data.HeteroData])
from torch_geometric.data import HeteroData
import torch_geometric.transforms as T
from torch_geometric.utils import convert
from torch_geometric.transforms import RandomLinkSplit
import torch.optim as optim
from torch_geometric.nn import ComplEx, DistMult, RotatE, TransE

"""
Goal: Test CNA feature extraction using KG node/edge embeddings

todo: build knowledge graph from patients metadata and the relationship of each sample with cna in genes

EFO extension by progenetix to describe cna events
id: EFO:0030063
label: copy number assessment
  |
  |-id: EFO:0030064
  | label: regional base ploidy
  |   |
  |   |-id: EFO:0030065
  |     label: copy-neutral loss of heterozygosity
  |
  |-id: EFO:0030066
    label: relative copy number variation
      |
      |-id: EFO:0030067
      | label: copy number loss
      |   |
      |   |-id: EFO:0030068
      |   | label: low-level copy number loss
      |   |
      |   |-id: EFO:0030069
      |     label: complete genomic deletion
      |
      |-id: EFO:0030070
        label: copy number gain
          |
          |-id: EFO:0030071
          | label: low-level copy number gain
          |
          |-id: EFO:0030072
             label: high-level copy number gain
             note: commonly but not consistently used for >=5 copies on a bi-allelic genome region
              |
              |-id: EFO:0030073
                 label: focal genome amplification
                 note: >-
                   commonly used for localized multi-copy genome amplification events where the
                   region does not extend >3Mb (varying 1-5Mb) and may exist in a large number of
                   copies
"""

class BuildCNAKg:
    def __init__(self, fout):
        self.opentarget_url = "https://api.platform.opentargets.org/api/v4/graphql"
        self.timeout = 30
        self.entities = {}
        self.relations = {}

        self.out = fout
        if( not os.path.exists(self.out) ):
            os.makedirs( self.out )

        self.processed = os.path.join( self.out, "kgdata" )
        if( not os.path.exists(self.processed) ):
            os.makedirs( self.processed )
    
    def test_fetch_opentarget_schema(self):
        # https://github.com/opentargets/open-targets-platform-mcp/blob/main/src/open_targets_platform_mcp/client/graphql.py

        # Create a client with schema fetching enabled
        async def async_conn():
            transport = AIOHTTPTransport( url = self.opentarget_url, headers = { "Content-Type": "application/json", }, timeout = self.timeout, )
            client = Client(transport=transport, fetch_schema_from_transport=True)
            async with client:
                # The schema is automatically fetched and stored in the client
                if not client.schema:
                    error_msg = "Failed to fetch schema from the GraphQL endpoint."
                    raise ValueError(error_msg)
                # Convert schema to SDL format
                return client.schema

        schema = asyncio.run( async_conn() )
        text = print_schema(schema)
        opath = os.path.join( self.out, "opentarget_schema.txt" )
        f = open(opath, 'w')
        f.write(text)
        f.close()

        #print_schema(client.schema)

    def __get_class_cnv(self, cnv):
        class_cnv = 'neutral'
        if( cnv >= 2.64 and cnv <= 3.36 ):
            class_cnv = 'gain'
        elif( cnv > 3.36 ):
            class_cnv = 'amp'
        elif( cnv >= 0.87 and cnv <= 1.32 ):
            class_cnv = 'loss'
        elif( cnv < 0.87 ):
            class_cnv = 'del'

        return class_cnv

    def _get_id(self, rtype, specific, value):
        if( rtype == "entity"):
            if(not specific in self.entities ):
                self.entities[specific] = {}
            if(not value in self.entities[specific] ):
                self.entities[specific][value] = len(self.entities[specific])
            return self.entities[specific][value]

        if( rtype == "relation"):
            if(not value in self.relations ):
                self.relations[value] = [ specific, len(self.relations) ]
            return self.relations[value]

    def get_patient_triplets_gdc(self, obj, splits):
        map_case = json.load( open( obj['map_case_clinfiles'], 'r') )
        patdata = json.load( open( obj['patient'], 'r') )
        for d in tqdm(patdata):
            uuid = d['uuid']
            case = map_case[uuid]
            _id = "Patient_%s" %( case )
            part = splits[case]
            
            triplets = set()
            kid = self._get_id('entity', 'patient', _id)
            for term in d:
                value = d[term]
                if(value not in ['', None, 'not_collected', 'not_provided']):
                    try:
                        value = value.lower()
                    except:
                        pass

                    if( (term.find('pathologic_stage')!=-1) or (term.find('cancer_stage')!=-1) ):
                        rel = "hasCancerStage"
                        rid = self._get_id('relation', 'patient|cancer_stage', rel)
                        oid = self._get_id('entity', 'cancer_stage', value)
                        triplets.add( "%s\t%s\t%s" %(_id, rel, value) )

                    elif( (term.find('gender')!=-1) or (term.find('sex')!=-1) ):
                        rel = "hasSex"
                        rid = self._get_id('relation', 'patient|sex', rel)
                        oid = self._get_id('entity', 'sex', value)
                        triplets.add( "%s\t%s\t%s" %(_id, rel, value) )

                    elif( (term.find('age')!=-1) and (term.find('stage')==-1) ):
                        rel = "hasAge"
                        rid = self._get_id('relation', 'patient|age', rel)
                        oid = self._get_id('entity', 'age', value)
                        triplets.add( "%s\t%s\t%s" %(_id, rel, value) )

                    elif( (term.find('ethnicity')!=-1) ):
                        rel = "hasEthnicity"
                        rid = self._get_id('relation', 'patient|ethnicity', rel)
                        oid = self._get_id('entity', 'ethnicity', value)
                        triplets.add( "%s\t%s\t%s" %(_id, rel, value) )

                    elif( (term.find('race')!=-1) ):
                        rel = "hasRace"
                        rid = self._get_id('relation', 'patient|race', rel)
                        oid = self._get_id('entity', 'race', value)
                        triplets.add( "%s\t%s\t%s" %(_id, rel, value) )

                    elif( (term.find('vital')!=-1) ):
                        rel = "hasVitalStatus"
                        rid = self._get_id('relation', 'patient|vital_status', rel)
                        oid = self._get_id('entity', 'vital_status', value)
                        triplets.add( "%s\t%s\t%s" %(_id, rel, value) )

                    elif( (term.find('tumor_tissue_site')!=-1) or (term.find('primary_site')!=-1) ):
                        rel = "hasTumorSite"
                        rid = self._get_id('relation', 'patient|tumor_site', rel)
                        oid = self._get_id('entity', 'tumor_site', value)
                        triplets.add( "%s\t%s\t%s" %(_id, rel, value) )
            
            lines = '\n'.join(triplets)
            opath = os.path.join(self.out, "kg_%s_%s.tsv" %(part, 'gdc') )
            opathall = os.path.join(self.out, "%s.tsv" %(part) )
            ls = [opath, opathall]
            for path in ls:
                f = open(path, 'a')
                f.write(lines+'\n')
                f.close()

    def get_patient_triplets_biostudies(self, obj, splits):
        patdata = obj['patient']
        for dpath in tqdm(os.listdir(patdata)):
            if( dpath.find("_metadata_")!=-1 and dpath.endswith('.json') ):
                study = dpath.replace('result_parsing_metadata_','').split('.')[0]
                dpath = os.path.join(patdata, dpath)
                d = json.load( open( dpath, 'r') )
                models = list(d['files'])
                part = splits[study]
                
                triplets = set()

                for uuid in models:
                    _id = "Patient_%s" %( study )
                    kid = self._get_id('entity', 'patient', _id)
                    for term in d["patient_tumor_info"]:
                        value = d["patient_tumor_info"][term]
                        if(value not in ['', None, 'not_collected', 'not_provided']):
                            try:
                                value = value.lower()
                            except:
                                pass
                            
                            if( (term.find('pathologic_stage')!=-1) or (term.find('cancer_stage')!=-1) ):
                                rel = "hasCancerStage"
                                rid = self._get_id('relation', 'patient|cancer_stage', rel)
                                oid = self._get_id('entity', 'cancer_stage', value)
                                triplets.add( "%s\t%s\t%s" %(_id, rel, value) )

                            elif( (term.find('gender')!=-1) or (term.find('sex')!=-1) ):
                                rel = "hasSex"
                                rid = self._get_id('relation', 'patient|sex', rel)
                                oid = self._get_id('entity', 'sex', value)
                                triplets.add( "%s\t%s\t%s" %(_id, rel, value) )

                            elif( (term.find('age')!=-1) and (term.find('stage')==-1) ):
                                rel = "hasAge"
                                rid = self._get_id('relation', 'patient|age', rel)
                                oid = self._get_id('entity', 'age', value)
                                triplets.add( "%s\t%s\t%s" %(_id, rel, value) )

                            elif( (term.find('ethnicity')!=-1) ):
                                rel = "hasEthnicity"
                                rid = self._get_id('relation', 'patient|ethnicity', rel)
                                oid = self._get_id('entity', 'ethnicity', value)
                                triplets.add( "%s\t%s\t%s" %(_id, rel, value) )

                            elif( (term.find('race')!=-1) ):
                                rel = "hasRace"
                                rid = self._get_id('relation', 'patient|race', rel)
                                oid = self._get_id('entity', 'race', value)
                                triplets.add( "%s\t%s\t%s" %(_id, rel, value) )

                            elif( (term.find('vital')!=-1) ):
                                rel = "hasVitalStatus"
                                rid = self._get_id('relation', 'patient|vital_status', rel)
                                oid = self._get_id('entity', 'vital_status', value)
                                triplets.add( "%s\t%s\t%s" %(_id, rel, value) )

                            elif( (term.find('tumor_tissue_site')!=-1) or (term.find('primary_site')!=-1) ):
                                rel = "hasTumorSite"
                                rid = self._get_id('relation', 'patient|tumor_site', rel)
                                oid = self._get_id('entity', 'tumor_site', value)
                                triplets.add( "%s\t%s\t%s" %(_id, rel, value) )
            
                lines = '\n'.join(triplets)
                opath = os.path.join(self.out, "kg_%s_%s.tsv" %(part, 'biostudies') )
                opathall = os.path.join(self.out, "%s.tsv" %(part) )
                ls = [opath, opathall]
                for path in ls:
                    f = open(path, 'a')
                    f.write(lines+'\n')
                    f.close()

    def _solve_cna_id(self, obj, source, fname):
        if( source == 'biostudies' ):
            _id = fname.replace('_processed.tsv','')
            case, uuid = _id.split('@')
        elif( source == 'gdc' ):
            map_case = json.load( open( obj['map_case_cnafiles'], 'r') )
            uuid = fname.replace('raw_','').split('.')[0]
            case = map_case[uuid]

        _id = "Patient_%s" %( case )
        return _id, case

    def get_cnaEvents_triplets(self, obj, splits, source):
        cnafolder = obj['cna']
        for fname in tqdm(os.listdir(cnafolder)):
            dpath = os.path.join(cnafolder, fname)
            patid, case = self._solve_cna_id(obj, source, fname)
            kid = self._get_id('entity', 'patient', patid)
            if(case in splits):
                part = splits[case]
                
                triplets = set()
                df = pd.read_csv(dpath, sep='\t')
                for i in df.index:
                    cnaid = str(uuid.uuid4())
                    _id = "cnaEvent_%s" %(cnaid)
                    kid = self._get_id('entity', 'cnv_event', _id)

                    gene = df.loc[i,'gene_name']
                    chrom = df.loc[i, 'chromosome']
                    cnv = df.loc[i, 'copy_number']
                    if( (str(cnv)!='nan') and (not str(cnv).startswith('2')) ):
                        rel = "belongsTo"
                        rid = self._get_id('relation', 'cnv_event|patient', rel)
                        triplets.add( "%s\t%s\t%s" %(_id, rel, patid) )

                        class_cnv = self.__get_class_cnv(cnv)
                        rel = "hasGene"
                        rid = self._get_id('relation', 'cnv_event|gene', rel)
                        oid = self._get_id('entity', 'gene', gene)
                        triplets.add( "%s\t%s\t%s" %(_id, rel, gene) )

                        rel = "hasChromosome"
                        rid = self._get_id('relation', 'cnv_event|chromosome', rel)
                        oid = self._get_id('entity', 'chromosome', chrom)
                        triplets.add( "%s\t%s\t%s" %(_id, rel, chrom) )

                        rel = "hasCnvNumber"
                        rid = self._get_id('relation', 'cnv_event|cnv_number', rel)
                        oid = self._get_id('entity', 'cnv_number', cnv)
                        triplets.add( "%s\t%s\t%.2f" %(_id, rel, cnv) )

                        rel = "hasCnvClass"
                        rid = self._get_id('relation', 'cnv_event|cnv_class', rel)
                        oid = self._get_id('entity', 'cnv_class', class_cnv)
                        triplets.add( "%s\t%s\t%s" %(_id, rel, class_cnv) )
                
                lines = '\n'.join(triplets)
                opath = os.path.join(self.out, "kg_%s_%s.tsv" %(part, source) )
                opathall = os.path.join(self.out, "%s.tsv" %(part) )
                ls = [opath, opathall]
                for path in ls:
                    f = open(path, 'a')
                    f.write(lines+'\n')
                    f.close()

        # test
        """
        sampleId123,hasSex,female
        sampleId123,hasAge,42
        sampleId123,hasCnaEvent,cnaevent1
        cnaevent1,hasGene,AC108479
        cnaevent1,hasScore,3
        cnaevent1,hasClassCna,gain

        """

    def __get_split_indexes(self, lst):
        lst = np.array(lst)
        np.random.shuffle(lst)
        idxs = np.array( list(range( len(lst) )) )
        tr, val, te = np.split(idxs, [int(.8 * len(idxs)), int(.9 * len(idxs))])
        splits = {}
        for i in tr:
            splits[ lst[i] ] = 'train'
        for i in val:
            splits[ lst[i] ] = 'valid'
        for i in te:
            splits[ lst[i] ] = 'test'
        return splits

    def get_partition_indexes_biostudy(self, obj):
        lst = []
        patdata = obj['patient']
        for dpath in os.listdir(patdata):
            if( dpath.find("_metadata_")!=-1 and dpath.endswith('.json') ):
                study = dpath.replace('result_parsing_metadata_','').split('.')[0]
                lst.append(study)

        splits = self.__get_split_indexes(lst)

        return splits

    def get_partition_indexes_gdc(self, obj):
        lst = set()
        map_case = json.load( open( obj['map_case_clinfiles'], 'r') )
        patdata = json.load( open( obj['patient'], 'r') )
        for d in patdata:
            uuid = d['uuid']
            case = map_case[uuid]
            lst.add(case)

        lst = list(lst)
        splits = self.__get_split_indexes(lst)

        return splits

    def export_relation_entities(self):
        opath = os.path.join(self.out, "entities.tsv")
        f = open(opath, 'w')
        f.write("node_type\tnumeric_id\tlabel\n")
        for typ in self.entities:
            for k in self.entities[typ]:
                v = self.entities[typ][k]
                f.write( "%s\t%i\t%s\n" %(typ, v, k) )
        f.close()

        opath = os.path.join(self.out, "relations.tsv")
        f = open(opath, 'w')
        f.write("relation_nodes\tnumeric_id\tlabel\n")
        for k in self.relations:
            els = self.relations[k]
            f.write( "%s\t%i\t%s\n" %(els[0], els[1], k) )
        f.close()

    def load_to_hetero_data(self):
        data = HeteroData()

        limit = 1000000
        ndc = { 'cnv_event': {} }
        ncnas = 0
        opath = os.path.join(self.out, "entities.tsv")
        i = 0
        f = open(opath, 'r')
        for line in f:
            if(i > 0):
                nt, nid, label = line.replace('\n','').split('\t')
                if( nt not in ndc ):
                    ndc[nt] = {}
                if( (nt != 'cnv_event') or len(ndc['cnv_event']) < limit ):
                    if( nt == 'cnv_number' ):
                        label = str( int( float(label) ) )
                    ndc[nt][label] = int(nid)

            i+=1
        f.close()
        
        '''
        https://github.com/davidlamprecht/AutoRDF2GML/blob/main/use-with-pyg/create-pyg-heterodata.py
        df = pd.read_csv(opath, sep='\t')
        for nt in df.node_type.unique():
            dfilt = df[ df['node_type'] == nt ]
            ndc[nt] = dict( zip( dfilt.label.values, dfilt.numeric_id.values ))
            data[nt] = torch.tensor( dfilt.numeric_id.values )
        '''
        
        opath = os.path.join(self.out, "relations.tsv")
        df = pd.read_csv(opath, sep='\t')
        rdc = dict( zip( df.label.values, df.relation_nodes.values ))
        ridc = dict( zip( df.label.values, df.numeric_id.values ))

        num_nodes = 0
        num_edges = 0
        edge_types = set()
        node_types = list(ndc)
        # feeding nodes
        for nt in ndc:
            node_ids = list(ndc[nt].values())
            data[nt].node_id = torch.tensor( node_ids )
            num_nodes += len(node_ids)
        
        # Feeding edges
        all_src = []
        all_rel_type = []
        all_dest = []
        rels = {}
        splits = ['train','valid','test']
        for part in splits:
            opathall = os.path.join(self.out, "%s.tsv" %(part) )
            f = open(opathall, 'r')
            for line in f:
                l = line.replace('\n','').split('\t')
                if( len(l) > 1 ):
                    s, p, o = l
                    nt1, nt2 = rdc[p].split('|')
                    pid = ridc[p]
                    
                    rid = "%s|%s|%s" %(nt1, p, nt2)
                    edge_types.add( (nt1, p, nt2) )

                    if( not rid in rels):
                        rels[rid] = { 'src': [], 'dest': [] }
                    
                    #print(nt1, nt2, s, o, type(o))
                    if( s in ndc[nt1] ):
                        #try:
                        rels[rid]['src'].append( ndc[nt1][s] )

                        if( nt2 == 'cnv_number' ):
                            o = str( int( float(o) ) )
                        rels[rid]['dest'].append( ndc[nt2][o] )
                        num_edges += 1

                        all_src.append(ndc[nt1][s])
                        all_rel_type.append(pid)
                        all_dest.append(ndc[nt2][o])

                        #except:
                        #    print(nt1, nt2, s, o)
            f.close()

        for rid in rels:
            nt1, p, nt2 = rid.split('|')
            e = ( '%s' %(nt1), '%s' %(p), '%s' %(nt2) )
            src = torch.tensor(rels[rid]['src'], dtype=torch.long)
            dest = torch.tensor(rels[rid]['dest'], dtype=torch.long)
            print( rid, rid.split('|'), len(src), len(dest) )
            data[ e ].edge_index = torch.stack( [ src, dest ], dim = 0 )

        data.node_types = node_types
        data.edge_types = list(edge_types)
        src = torch.tensor( all_src, dtype=torch.long)
        dest = torch.tensor( all_dest, dtype=torch.long)
        data.edge_index = torch.stack( [ src, dest ], dim = 0 )
        data.edge_type = torch.tensor(all_rel_type)

        data = T.ToUndirected()(data)
        opath = os.path.join(self.processed, "cnakg")
        torch.save(data, opath)

    def build_kg(self, dts):
        '''
        for dt in dts:
            obj = dts[dt]
            source = obj['sourcedb']
            if(source == 'gdc'):
                print('Extracting from gdc')
                splits = self.get_partition_indexes_gdc(obj)
                self.get_patient_triplets_gdc(obj, splits)
                self.get_cnaEvents_triplets(obj, splits, source)

            if(source == 'biostudies'):
                print('Extracting from biostudies')
                splits = self.get_partition_indexes_biostudy(obj)
                self.get_patient_triplets_biostudies(obj, splits)
                self.get_cnaEvents_triplets(obj, splits, source)

        print('Exporting relation')
        self.export_relation_entities()
        '''
        
        print('Loading to torch hetero data')
        self.load_to_hetero_data()

        """
        Statistics:
        patient|hasRace|race ['patient', 'hasRace', 'race'] 285 285
        patient|hasAge|age ['patient', 'hasAge', 'age'] 519 519
        patient|hasCancerStage|cancer_stage ['patient', 'hasCancerStage', 'cancer_stage'] 504 504
        patient|hasSex|sex ['patient', 'hasSex', 'sex'] 584 584
        patient|hasVitalStatus|vital_status ['patient', 'hasVitalStatus', 'vital_status'] 458 458
        patient|hasEthnicity|ethnicity ['patient', 'hasEthnicity', 'ethnicity'] 329 329
        patient|hasTumorSite|tumor_site ['patient', 'hasTumorSite', 'tumor_site'] 593 593
        cnv_event|hasCnvNumber|cnv_number ['cnv_event', 'hasCnvNumber', 'cnv_number'] 293307 293307
        cnv_event|hasGene|gene ['cnv_event', 'hasGene', 'gene'] 293307 293307
        cnv_event|belongsTo|patient ['cnv_event', 'belongsTo', 'patient'] 293307 293307
        cnv_event|hasChromosome|chromosome ['cnv_event', 'hasChromosome', 'chromosome'] 293307 293307
        cnv_event|hasCnvClass|cnv_class ['cnv_event', 'hasCnvClass', 'cnv_class'] 293307 293307
        """


    def _train(self, model, loader, optimizer):
        model.train()
        total_loss = total_examples = 0
        for head_index, rel_type, tail_index in loader:
            optimizer.zero_grad()
            loss = model.loss(head_index, rel_type, tail_index)
            loss.backward()
            optimizer.step()
            total_loss += float(loss) * head_index.numel()
            total_examples += head_index.numel()
        return total_loss / total_examples


    @torch.no_grad()
    def _test(self, model, data):
        model.eval()
        return model.test(
            head_index=data.edge_index[0],
            rel_type=data.edge_type,
            tail_index=data.edge_index[1],
            batch_size=20000,
            k=10,
        )
    
    def _reduce_dataset(self, data, keyfilter='cnv_event', nsamples=50000):
        opath = os.path.join(self.out, "relations.tsv")
        df = pd.read_csv(opath, sep='\t')
        rdc = dict( zip( df.label.values, df.relation_nodes.values ))
        ridc = dict( zip( df.label.values, df.numeric_id.values ))
        
        all_src = []
        all_rel_type = []
        all_dest = []

        tmp = data
        for et in data.edge_types:
            s, p, o = et
            if( ((s == keyfilter) or (p == keyfilter) or (o == keyfilter)) ):
                src = data[et].edge_index[0]
                dest = data[et].edge_index[1]
                total = len(src)
                if(nsamples < total ):
                    idxs = np.random.choice( total, nsamples, replace=False)
                    tmp[et].edge_index = torch.stack( [ src[idxs], dest[idxs] ], dim = 0 )
            
            all_src.extend( tmp[et].edge_index[0].tolist() )
            all_rel_type.extend( [ ridc[p.replace('rev_','')] ]*nsamples )
            all_dest.extend( tmp[et].edge_index[1].tolist() )

        all_src = torch.tensor(all_src)
        all_rel_type = torch.tensor(all_rel_type)
        all_dest = torch.tensor(all_dest)
        src = torch.tensor( all_src, dtype=torch.long)
        dest = torch.tensor( all_dest, dtype=torch.long)
        tmp.edge_index = torch.stack( [ src, dest ], dim = 0 )
        tmp.edge_type = torch.tensor(all_rel_type)

        return tmp

    def train_kg(self):
        """
        https://github.com/pyg-team/pytorch_geometric/blob/master/examples/kge_fb15k_237.py
        
        Mean Reciprocal Rank (MRR) and Mean Rank (MR) are the most common rank-based metrics used to evaluate Knowledge Graph Embedding (KGE) models, particularly for link prediction tasks (e.g., predicting missing entities or relations in a triple (h, r, t). They assess how high a model ranks the correct ("ground-truth") entity compared to all other candidates. 
        Key Metrics Explained
        Mean Reciprocal Rank (MRR): Measures the average of the inverse ranks of the first relevant item across all queries.
        Formula: 
        .
        Interpretation: Ranging from 0 to 1, a higher MRR (closer to 1) indicates that the correct entity is consistently ranked near the top of the list.
        Preference: MRR is often preferred over Mean Rank as it is more robust to outliers.
        Mean Rank (MR): Calculates the average position of the ground-truth entity.
        Interpretation: A lower MR (closer to 1) is better.
        Drawback: MR is highly sensitive to outliers, meaning a single, extremely poor prediction can heavily skew the average.
        Hits@K (Hits at K): Measures the proportion of correct entities ranked in the top positions (commonly K=1,2,3...10)

        """
        all_edges = [ ('patient', 'hasRace', 'race'), ('patient', 'hasAge', 'age'), ('patient', 'hasCancerStage', 'cancer_stage'), ('patient', 'hasSex', 'sex'), ('patient', 'hasVitalStatus', 'vital_status'), ('patient', 'hasEthnicity', 'ethnicity'), ('patient', 'hasTumorSite', 'tumor_site'), ('cnv_event', 'hasCnvNumber', 'cnv_number'), ('cnv_event', 'hasGene', 'gene'), ('cnv_event', 'belongsTo', 'patient'), ('cnv_event', 'hasChromosome', 'chromosome'), ('cnv_event', 'hasCnvClass', 'cnv_class') ]
        rev_all_edges = [ ('race', 'rev_hasRace', 'patient'), ('age', 'rev_hasAge', 'patient'), ('cancer_stage', 'rev_hasCancerStage', 'patient'), ('sex', 'rev_hasSex', 'patient'), ('vital_status', 'rev_hasVitalStatus', 'patient'), ('ethnicity', 'rev_hasEthnicity', 'patient'), ('tumor_site', 'rev_hasTumorSite', 'patient'), ('cnv_number', 'rev_hasCnvNumber', 'cnv_event'), ('gene', 'rev_hasGene', 'cnv_event'), ('patient', 'rev_belongsTo', 'cnv_event'), ('chromosome', 'rev_hasChromosome', 'cnv_event'), ('cnv_class', 'rev_hasCnvClass', 'cnv_event') ]
        sel_edges = [ ('patient', 'hasTumorSite', 'tumor_site'), ('cnv_event', 'hasCnvNumber', 'cnv_number'), ('cnv_event', 'hasGene', 'gene'), ('cnv_event', 'belongsTo', 'patient'), ('cnv_event', 'hasCnvClass', 'cnv_class') ]
        rev_sel_edges = [ ('tumor_site', 'rev_hasTumorSite', 'patient'), ('cnv_number', 'rev_hasCnvNumber', 'cnv_event'), ('gene', 'rev_hasGene', 'cnv_event'), ('patient', 'rev_belongsTo', 'cnv_event'), ('cnv_class', 'rev_hasCnvClass', 'cnv_event') ]

        # Assume 'data' is your HeteroData object
        models = list(model_map)
        models = ['rotate']
        nepochs = 20
        key_filter = 'cnv_event'
        nreduction = 20000
        
        transform = RandomLinkSplit(
            num_val = 0.1,  # 10% validation
            num_test = 0.1, # 10% test
            edge_types = sel_edges, # Specific edge to split
            rev_edge_types = rev_sel_edges, # Reverse edge
            is_undirected = False,
            add_negative_train_samples = True,
            split_labels = True
        )

        print('Loading data')
        path = os.path.join(self.processed, "cnakg")
        data = torch.load(path, weights_only=False)
        
        print("Reducing data")
        subdata = self._reduce_dataset( data, keyfilter = key_filter, nsamples = nreduction)

        print('Transforming and spliting')
        train_data, val_data, test_data = transform(subdata)

        model_map = {
            'transe': TransE,
            'complex': ComplEx,
            'distmult': DistMult,
            'rotate': RotatE,
        }
        model_arg_map = {'rotate': {'margin': 9.0}}

        opath = os.path.join(self.out, "report.txt")
        f = open(opath, 'w')
        for m in models:
            line = f'In {m}\n'
            print( line )
            f.write(line)

            num_edge_types = len(data.edge_types)
            model = model_map[m](
                num_nodes = train_data.num_nodes,
                num_relations = train_data.num_edges,
                hidden_channels = 50,
                **model_arg_map.get(m, {}),
            ).to(device)

            optimizer_map = {
                'transe': optim.Adam(model.parameters(), lr=0.01),
                'complex': optim.Adagrad(model.parameters(), lr=0.001, weight_decay=1e-6),
                'distmult': optim.Adam(model.parameters(), lr=0.0001, weight_decay=1e-6),
                'rotate': optim.Adam(model.parameters(), lr=1e-3),
            }

            optimizer = optimizer_map[m]

            before = datetime.datetime.now()

            loader = model.loader(
                head_index = train_data.edge_index[0],
                rel_type = train_data.edge_type,
                tail_index = train_data.edge_index[1],
                batch_size = 1000,
                shuffle = True,
            )

            for epoch in range(1, nepochs):
                loss = self._train(model, loader, optimizer)
                print(f'\tEpoch: {epoch:03d}, Loss: {loss:.4f}')
                if epoch % 25 == 0:
                    rank, mrr, hits = self._test(model, val_data)
                    line = f'\tEpoch: {epoch:03d}, Val Mean Rank: {rank:.2f}, Val MRR: {mrr:.4f}, Val Hits@10: {hits:.4f}\n'
                    print(line)
                    f.write(line)


            after = datetime.datetime.now()  
            
            datefmt = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            delta_time = after-before
            formatted = str(datetime.timedelta(seconds = delta_time.seconds))

            rank, mrr, hits_at_10 = self._test(test_data)
            line = f'\tTest Mean Rank: {rank:.2f}, Test MRR: {mrr:.4f}, Test Hits@10: {hits_at_10:.4f}\n'
            line += "\tExecution time: %s\n" %(formatted)
            print(line)
            f.write(line)

            opath = os.path.join(self.processed, "model_%s.pt" %(m) )
            torch.save(model, opath)
        f.close()


    def call_train(self):
        models = ['TransE', 'TransE_l1', 'TransE_l2', 'TransR', 'RESCAL', 'DistMult', 'ComplEx', 'RotatE', 'SimplE']
        for m in models:
            fout = os.path.join(self.out, m)
            if(not os.path.exists(fout)):
                os.makedirs(fout)

            cmd = "dglke_train --model_name %s --data_path %s --format udd_hrt --save_path %s --delimiter '\t' --batch_size 400 --neg_sample_size 200 --hidden_dim 400 --gamma 19.9 --lr 0.25 --max_step 3000 --log_interval 100 --batch_size_eval 16 --test -adv --regularization_coef 1.00E-09 --num_thread 1 --num_proc 4" %( m, self.out, fout )

            cmd = cmd.split(" ")
            out = run( cmd, stdout=PIPE, stderr=PIPE )
            if( out.returncode != 0 ):
                raise Exception(out.stderr.decode())
            else:
                pass

    def get_gene_association_scores(self):
        # The efo id for colorectal cancer is EFO_0004142
        disease = "EFO_0004142"
        score_sources = {'gwas_credible_sets': 1, 'gene_burden': 1, 'eva': 1, 'genomics_england': 1, 'gene2phenotype': 1, 'uniprot_literature': 1, 'uniprot_variants': 1, 'orphanet': 1, 'clingen': 1, 'cancer_gene_census': 1, 'intogen': 1, 'eva_somatic': 1, 'cancer_biomarkers': 1, 'chembl': 1, 'crispr_screen': 1, 'crispr': 1, 'reactome': 1, 'europepmc': 0.2, 'expression_atlas': 0.2, 'impc': 0.2, 'ot_crispr_validation': 0.5, 'ot_crispr': 0.5, 'encore': 0.5}
        datasources = []
        for k in score_sources:
            datasources.append( { "id": k, "weight": score_sources[k], "propagate": True, "required": False } )

        variables = {
          "id": disease,
          "index": 0,
          "size": 500,
          "sortBy": "score",
          "enableIndirect": True, # retrieve the top
          "entity": "disease",
          "entitySearch": "",
          "datasources": datasources
        }

        query = """
query DiseaseAssociationsQuery(
  $id: String!
  $index: Int!
  $size: Int!
  $sortBy: String!
  $enableIndirect: Boolean!
  $datasources: [DatasourceSettingsInput!]
  $rowsFilter: [String!]
  $facetFilters: [String!]
  $entitySearch: String!
) {
  disease(efoId: $id) {
    id
    name
    associatedTargets(
      page: { index: $index, size: $size }
      orderByScore: $sortBy
      enableIndirect: $enableIndirect
      datasources: $datasources
      Bs: $rowsFilter
      facetFilters: $facetFilters
      BFilter: $entitySearch
    ) {
      count
      rows {
        target {
          id
          approvedSymbol
          approvedName
          prioritisation {
            items {
              key
              value
            }
          }
        }
        score
        datasourceScores {
          componentId: id
          score
        }
      }
    }
  }
}

        """
        query = gql(query)
        # Create a client with schema fetching enabled
        async def async_conn():
            transport = AIOHTTPTransport( url = self.opentarget_url, headers = { "Content-Type": "application/json", }, timeout = self.timeout, )
            client = Client(transport=transport, fetch_schema_from_transport=True)
            result = await client.execute_async(query, variable_values=variables)

            return result

        result = asyncio.run( async_conn() )
        opath = os.path.join( self.processed, "top_gene_assoc.json" )
        json.dump(result, open(opath, 'w') )

        variables["enableIndirect"] = False
        result = asyncio.run( async_conn() )
        opath = os.path.join( self.processed, "bottom_gene_assoc.json" )
        json.dump(result, open(opath, 'w') )

    
    def _get_node_map(self, ntypes):
        path = os.path.join(self.out, "entities.tsv")
        ndc = {}
        i = 0
        f = open(path, 'r')
        for line in f:
            if(i > 0):
                nt, nid, label = line.replace('\n','').split('\t')
                if( nt in ntypes ):
                    _id = "%s_%s" %(nt, nid)
                    ndc[_id] = label
            i+=1
        f.close()

        return ndc

    def _prepare_networkx_graph(self, data):
        g = nx.Graph()
        
        opath = os.path.join(self.processed, "chakg_networkx.graphml")
        if( not os.path.exists(opath) ):
            node_types = ['patient', 'gene', 'race', 'age', 'sex', 'ethnicity', 'vital_status', 'tumor_site', 'cancer_stage', 'cnv_class', 'cnv_number']
            mp = self._get_node_map(node_types)
            edges = [ 
                ('patient', 'hasRace', 'race'), 
                ('patient', 'hasAge', 'age'), 
                ('patient', 'hasCancerStage', 'cancer_stage'), 
                ('patient', 'hasSex', 'sex'), 
                ('patient', 'hasVitalStatus', 'vital_status'), 
                ('patient', 'hasEthnicity', 'ethnicity'), 
                ('patient', 'hasTumorSite', 'tumor_site'), 
                #('cnv_event', 'hasCnvNumber', 'cnv_number'), 
                ('cnv_event', 'hasGene', 'gene'), 
                ('cnv_event', 'belongsTo', 'patient'), 
                #('cnv_event', 'hasChromosome', 'chromosome'), 
                ('cnv_event', 'hasCnvClass', 'cnv_class')
            ]

            cnv_numbers = list( map( lambda x: 'cnv_number_%i' %(x), data[ ('cnv_event', 'hasCnvNumber', 'cnv_number') ].edge_index[1].tolist() ))
            
            for et in tqdm(edges):
                print(et)
                s, p, o = et
                for i, v in enumerate( data[et].edge_index[0] ):
                    src_val = data[ et ].edge_index[0][i]
                    src = "%s_%i" %(s, src_val)
                    if(s != 'cnv_event'):
                        src = mp[ src ]

                    dest_val = data[ et ].edge_index[1][i]
                    dest = "%s_%i" %(o, dest_val)
                    if(o != 'cnv_event'):
                        dest = mp[ dest ]
                    
                    if(s == 'cnv_event'):
                        cnv_number = float( mp[ cnv_numbers[i] ] )
                        g.add_edge(src, dest, weight = cnv_number, edge_type=p)
                    else:
                        g.add_edge(src, dest, edge_type=p)
                    g.nodes()[src]['type'] = s
                    g.nodes()[dest]['type'] = o
            
            nx.write_graphml(g, opath)
        else:
            g = nx.read_graphml(opath)
        
        """
        patient = list( map( lambda x: 'patient_%i' %(x), data[ ('cnv_event', 'belongsTo', 'patient') ].edge_index[1].tolist() ))
        for i, c in enumerate(patient):
            pid = patient[i]

            gene = data[ ('cnv_event', 'hasGene', 'gene') ].edge_index[1][i]
            gene = mp[ "gene_%i" %(gene) ]

            class_cnv = data[ ('cnv_event', 'hasCnvClass', 'cnv_class') ].edge_index[1][i]
            class_cnv = mp[ "cnv_class_%i" %(class_cnv) ]
            
            cnv_number = data[ ('cnv_event', 'hasCnvNumber', 'cnv_number') ].edge_index[1][i]
            cnv_number = float(mp[ "cnv_number_%i" %(cnv_number) ])
            
            g.add_edge(pid, gene, weight = cnv_number, edge_type='hasGene')
            g.nodes()[gene]['type'] = 'gene'
            g.add_edge(pid, class_cnv, weight = cnv_number, edge_type='hasCnvClass')
            g.nodes()[class_cnv]['type'] = 'cnv_class'

            g.nodes()[pid]['type'] = 'patient'
        """


        return g

    def prepare_metapath_input(self, y_info_type = 'opentarget_assoc', y_data = {}):
        print('Loading data')
        path = os.path.join(self.processed, "cnakg")
        data = torch.load(path, weights_only=False)
        
        metapath = [
            ('gene', 'rev_hasGene', 'cnv_event'), 
            ('cnv_event', 'belongsTo', 'patient'), 
            ('patient', 'hasAge', 'age')
        ]

        edges = [ 
            ('patient', 'hasRace', 'race'), 
            ('patient', 'hasAge', 'age'), 
            ('patient', 'hasCancerStage', 'cancer_stage'), 
            ('patient', 'hasSex', 'sex'), 
            ('patient', 'hasVitalStatus', 'vital_status'), 
            ('patient', 'hasEthnicity', 'ethnicity'), 
            ('patient', 'hasTumorSite', 'tumor_site'), 
            #('cnv_event', 'hasCnvNumber', 'cnv_number'), 
            ('cnv_event', 'hasGene', 'gene'), 
            ('gene', 'rev_hasGene', 'cnv_event'), 
            ('cnv_event', 'belongsTo', 'patient'), 
            #('cnv_event', 'hasChromosome', 'chromosome'), 
            ('cnv_event', 'hasCnvClass', 'cnv_class')
        ]
        edge_dict = {}
        for e in edges:
            edge_dict[e] = data[e]['edge_index']

        opath = os.path.join(self.processed, y_info_type+"_cnakg_gene_labels")
        if( not os.path.exists(opath) ):
            labels = self._get_node_map( ['gene'] )
            #y_index = data['gene']['node_id']
            y_index = []
            
            if( y_info_type == 'opentarget_assoc' ):
                y_info = { 'y_ot_assoc': [] }
                for i in data['gene']['node_id'].tolist():
                    gn = labels['gene_'+str(i)]
                    if( (gn in y_data) ):
                        v = y_data[gn]
                        y_info['y_ot_assoc'].append(v)
                        y_index.append(i)

            elif( y_info_type == 'civic' ):
                # Adjust gene node labels to civic information
                y_info = {}
                evtypes = set()
                sigs = set()
                for g in y_data:
                    evtypes.update( list( y_data[g]['evtype']) )
                    sigs.update( list( y_data[g]['significance']) )

                for i in data['gene']['node_id'].tolist():
                    gn = labels['gene_'+str(i)]
                    for e in evtypes:
                        k = "y_"+e
                        if( not k in y_info ):
                            y_info[k] = []
                        v = 0
                        if( (gn in y_data) and (e in y_data[gn]['evtype']) ):
                            v = 1
                            y_info[k].append(v)
                            y_index.append(i)

                    for e in sigs:
                        k = "y_"+e
                        if( not k in y_info ):
                            y_info[k] = []
                        v = 0
                        if( (gn in y_data) and (e in y_data[gn]['significance']) ):
                            v = 1
                            y_info[k].append(v)
                            y_index.append(i)

            data['gene']['y_index'] = torch.tensor(y_index, dtype=torch.long)
            for k in y_info:
                data['gene'][k] = torch.tensor(y_info[k], dtype=torch.long)
            torch.save(data, opath)
        else:
            data = torch.load(opath, weights_only=False)

        return metapath, edge_dict, data

    def  extract_gene_networkx_topology_metrics(self):
        # Apart from degree, all the other metrics takes a lot of time to compute in a normal pc

        print('Loading data')
        path = os.path.join(self.processed, "cnakg")
        data = torch.load(path, weights_only=False)

        # = convert.to_networkx(data)
        
        g = self._prepare_networkx_graph(data)

        gids = g.nodes
        gdc = g.nodes()
        gene_nodes = list( filter( lambda n: gdc[n]['type'] == 'gene', gids ) )
        
        measures = ['degree', 'eigenvector', 'katz', 'closeness', 'information', 'betweenness', 'clustering']
        measures = ['degree', 'abs_degree']
        res = { 'degree': nx.degree_centrality(g), 'abs_degree': nx.degree(g) }
        features = {}
        for n in tqdm(gene_nodes):
            features[n] = { }
            for m in measures:
                features[n][m] = res[m][n]
        """
        features = {}
        for n in tqdm(gene_nodes):
            features[n] = { }
            for m in measures:
                if(m != 'clustering'):
                    features[n][m] = eval("nx.%s_centrality(g)[n]" %(m))
                    if( m == 'degree' ):
                        features[n]['abs_degree'] = eval("nx.%s(g)[n]" %(m))
                else:
                    features[n][m] = nx.clustering(g, n)
        """

        opath = os.path.join(self.processed, 'topology_features.tsv')
        header = ['gene'] + measures
        lines = [ header ]
        for g in features:
            el = [g] + list(features[g].values())
            lines.append(el)
        lines = list( map( lambda x: '\t'.join( [ str(y) for y in x ] ), lines ))
        f = open( opath, "w")
        f.write("\n".join(lines) + "\n")
        f.close()

    def _extract_features_opentarget_sources(self):
        opath = os.path.join( self.processed, "top_gene_assoc.json" )
        d = json.load( open(opath, 'r') )
        rs = d['disease']['associatedTargets']['rows']
        opath = os.path.join( self.processed, "bottom_gene_assoc.json" )
        d = json.load( open(opath, 'r') )
        rs += d['disease']['associatedTargets']['rows']

        features = set()
        dt = {}
        for r in rs:
            keys = ['y'] + list( map( lambda x: x['componentId'], r['datasourceScores'] ))
            features.update(keys)
            values = [ r['score'] ] + list( map( lambda x: x['score'], r['datasourceScores'] ))
            dc = dict( zip( keys, values ) )
            dt[ r['target']['approvedSymbol'] ] = dc

        opath = os.path.join(self.processed, 'opentarget_features.tsv')
        header = ['gene'] + list(features)
        lines = [ header ]
        for g in dt:
            values = []
            for f in features:
                try:
                    v = dt[g][f]
                except:
                    v = 0
                values.append(v)

            el = [g] + values
            lines.append(el)
        lines = list( map( lambda x: '\t'.join( [ str(y) for y in x ] ), lines ))
        f = open( opath, "w")
        f.write("\n".join(lines) + "\n")
        f.close()

    def merge_topology_opentarget_featSelection_regression(self):
        # Reference: https://www.geeksforgeeks.org/machine-learning/how-to-perform-feature-selection-for-regression-data/
        
        path = os.path.join(self.processed, 'topology_features.tsv')
        tmp = pd.read_csv(path, sep='\t')
        genes = set(tmp.gene.values)
        tdc = {}
        for i in tmp.index:
            gene = tmp.loc[i, 'gene']
            dg = tmp.loc[i, 'degree']
            adg = tmp.loc[i, 'abs_degree']
            tdc[gene] = [dg, adg]

        path = os.path.join(self.processed, 'opentarget_features.tsv')
        df = pd.read_csv(path, sep='\t')
        df = df[ df.gene.isin(genes) ]
        print('Remaining after opentarget filtering:', len(df))
        x_cols = list( filter( lambda x: (x not in [ 'y', 'gene']), list(df.columns) ))
        df['degree'] = [ tdc[g][0] for g in df.gene ]
        df['abs_degree'] = [ tdc[g][1] for g in df.gene ]
        X = df[ x_cols ]#.values
        y = df['y'].values

        from sklearn.model_selection import train_test_split
        from sklearn.linear_model import LinearRegression
        from sklearn.feature_selection import RFE
        from sklearn.metrics import mean_squared_error

        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
        model = LinearRegression()
        rfe = RFE(model, n_features_to_select=5)
        rfe = rfe.fit(X_train, y_train)
        selected_features = X_train.columns[rfe.support_]
        print("Selected Features:", selected_features)

        X_train_rfe = rfe.transform(X_train)
        X_test_rfe = rfe.transform(X_test)
        model.fit(X_train_rfe, y_train) 

        y_pred = model.predict(X_test_rfe)
        mse = mean_squared_error(y_test, y_pred)
        print("Mean Squared Error:", mse)

        """
        Remaining after opentarget filtering: 494
        Selected Features: Index(['clingen', 'intogen', 'eva', 'chembl', 'uniprot_literature'], dtype='object')
        Mean Squared Error: 0.004246576898605301

        Remaining after opentarget filtering: 871
        Selected Features: Index(['uniprot_variants', 'chembl', 'gwas_credible_sets', 'clingen',
               'cancer_gene_census'],
              dtype='object')
        Mean Squared Error: 0.005905122226428685

        """

    def run(self):
        #self.test_fetch_opentarget_schema()
        datasources = {
            "tcga_coad": {
                "sourcedb": "gdc",
                "cna": "/var/www/html/gdc_explorer_app/data_processed/TCGA-COAD/copy-number-variation/files/",
                "patient": "/var/www/html/gdc_explorer_app/data_processed/TCGA-COAD/clinical/data_cases.json",
                "map_case_clinfiles": "/var/www/html/gdc_explorer_app/data_processed/TCGA-COAD/clinical/mapping_file_case.json",
                "map_case_cnafiles": "/var/www/html/gdc_explorer_app/data_processed/TCGA-COAD/copy-number-variation/mapping_file_case.json"
            },
            "biostudies": {
                "sourcedb": "biostudies",
                "cna": "/mnt/yasdata/home/yasmmin/Dropbox/portfolio_2025/gdc_explorer_app_web/gdc_exploration_project/copy_number_variation/biostudies/processed/",
                "patient": "/mnt/yasdata/home/yasmmin/Dropbox/portfolio_2025/gdc_explorer_app_web/gdc_exploration_project/copy_number_variation/biostudies/"
            }

        }
        #self.build_kg(datasources)
        #self.train_kg()
        
        #self.get_gene_association_scores()
        #self.extract_gene_networkx_topology_metrics()
        
        #self._extract_features_opentarget_sources()
        self.merge_topology_opentarget_featSelection_regression()

        #self.call_train()
        # dglke is no longer available and it is not in line with the modern methods - Substitute by torch geometric
        # model - https://github.com/pyg-team/pytorch_geometric/blob/master/examples/kge_fb15k_237.py
        # Create torch geom custom Dataset - https://pytorch-geometric.readthedocs.io/en/latest/tutorial/load_csv.html

if( __name__ == "__main__" ):
    out = '/mnt/yasdata/home/yasmmin/Dropbox/portfolio_2025/gdc_explorer_app_web/gdc_exploration_project/copy_number_variation/cnakg'
    
    o = BuildCNAKg( out )
    o.run()