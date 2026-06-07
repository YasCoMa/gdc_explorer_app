// Test mljs

function shuffleArray(array) {
    for (var i = array.length - 1; i > 0; i--) {
        var j = Math.floor(Math.random() * (i + 1));
        var temp = array[i];
        array[i] = array[j];
        array[j] = temp;
    }
}

function error(predicted, expected) {
    let misclassifications = 0;
    for (var index = 0; index < predicted.length; index++) {
        console.log(`truth : ${expected[index]} and prediction : ${predicted[index]}`);
        if (predicted[index] !== expected[index]) {
            misclassifications++;
        }
    }
    return misclassifications;
}

function test() {
    const result = clf.predict(testSetX);
    const testSetLength = testSetX.length
    const predictionError = error(result, testSetY);
    console.log(`Test Set Size = ${testSetLength} and number of Misclassifications = ${predictionError}`);
}

r = await fetch('https://raw.githubusercontent.com/mljs/ml/refs/heads/main/examples/leafDataset/leaf.csv')
t = await r.text()

let types = new Set(); // To gather UNIQUE classes
df.forEach((row) => {
    types.add(row[0]);
});
typesArray = [...types]; // To save the different types of classes.

df = t.split('\n').slice(1,-1).map( e => e.split(',') )
shuffleArray(df)

X = df.map( e => e.slice(2, 16).map( e => parseFloat(e) ) )
y = df.map( e => typesArray.indexOf(e[0]) )

seperationSize = 0.9 * df.length;
trainingSetX = X.slice(0, seperationSize);
trainingSetY = y.slice(0, seperationSize);
testSetX = X.slice(seperationSize);
testSetY = y.slice(seperationSize);

clf = new ML.NaiveBayes.GaussianNB();
clf.train(trainingSetX, trainingSetY)
results = nb.predict(testSetX)

