import os.path as osp

import json
import torch
import torch_geometric
device = 'cuda' if torch.cuda.is_available() else 'cpu'
torch.serialization.add_safe_globals([torch_geometric.data.hetero_data.HeteroData])

from torch_geometric.datasets import AMiner
from torch_geometric.nn import MetaPath2Vec

from test_kggraph_embedding import BuildCNAKg

class EmbeddingMetapath:
    def __init__(self, fout):
        self.bkg = BuildCNAKg(fout)

    def get_inputs(self):
        path = '/var/www/html/gdc_explorer_app/data_processed/TCGA-COAD/copy-number-variation/y_civic.json'
        dc = json.load( open(path, 'r') )
        #edges, edge_dict, data = self.bkg.prepare_metapath_input('civic', dc)

        path = '/mnt/yasdata/home/yasmmin/Dropbox/portfolio_2025/gdc_explorer_app_web/gdc_exploration_project/copy_number_variation/cnakg/kgdata/top_gene_assoc.json'
        d = json.load( open(path, 'r') )
        rs = d['disease']['associatedTargets']['rows']
        
        path = '/mnt/yasdata/home/yasmmin/Dropbox/portfolio_2025/gdc_explorer_app_web/gdc_exploration_project/copy_number_variation/cnakg/kgdata/bottom_gene_assoc.json'
        d = json.load( open(path, 'r') )
        rs += d['disease']['associatedTargets']['rows']
        dc = {}
        for r in rs:
            gene = r['target']['approvedSymbol']
            score = r['score']
            dc[gene] = 1
            if(score < 0.4):
                dc[gene] = 0

        metapath, edge_dict, data = self.bkg.prepare_metapath_input('opentarget_assoc', dc)
        
        return metapath, edge_dict, data

    def setup_model(self, metapath, edge_index_dict, embdim = 128):
        wlen = len(metapath)
        ctxsize = wlen+1
        model = MetaPath2Vec( edge_index_dict, embedding_dim=embdim,
                     metapath=metapath, walk_length=wlen, context_size=ctxsize,
                     walks_per_node=2, num_negative_samples=5,
                     sparse=True).to(device)

        loader = model.loader(batch_size=128, shuffle=True, num_workers=6)
        optimizer = torch.optim.SparseAdam(list(model.parameters()), lr=0.01)

        return model, loader, optimizer

    def train(self, epoch, model, loader, optimizer, data, node, y_layer, log_steps=100, eval_steps=2000):
        model.train()

        total_loss = 0
        for i, (pos_rw, neg_rw) in enumerate(loader):
            optimizer.zero_grad()
            loss = model.loss(pos_rw.to(device), neg_rw.to(device))
            loss.backward()
            optimizer.step()

            total_loss += loss.item()
            if (i + 1) % log_steps == 0:
                print(f'Epoch: {epoch}, Step: {i + 1:05d}/{len(loader)}, '
                      f'Loss: {total_loss / log_steps:.4f}')
                total_loss = 0

            if (i + 1) % eval_steps == 0:
                acc = self.test(model, data, 'gene', 'y_ot_assoc')
                print(f'Epoch: {epoch}, Step: {i + 1:05d}/{len(loader)}, '
                      f'Acc: {acc:.4f}')


    @torch.no_grad()
    def test( self, model, data, node, y_layer, train_ratio=0.1):
        model.eval()

        bat = None
        bat = data[node]['y_index']
        z = model(node, batch = bat )
        y = data[node][y_layer]
        print('z', len(z), 'y', len(y) ) #it is finally working
        
        perm = torch.randperm(z.size(0))
        train_perm = perm[:int(z.size(0) * train_ratio)]
        test_perm = perm[int(z.size(0) * train_ratio):]

        return model.test(z[train_perm], y[train_perm], z[test_perm], y[test_perm],
                          max_iter=150)

    def run(self):
        metapath, edge_index_dict, data = self.get_inputs()
        model, loader, optimizer = self.setup_model(metapath, edge_index_dict)

        for epoch in range(1, 200):
            self.train(epoch, model, loader, optimizer, data, 'gene', 'y_ot_assoc')
            acc = self.test(model, data, 'gene', 'y_ot_assoc')
            print(f'Epoch: {epoch}, Accuracy: {acc:.4f}')

if( __name__ == "__main__" ):
    out = '/mnt/yasdata/home/yasmmin/Dropbox/portfolio_2025/gdc_explorer_app_web/gdc_exploration_project/copy_number_variation/cnakg'
    
    o = EmbeddingMetapath( out )
    o.run()