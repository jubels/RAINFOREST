# RAINFOREST
RAINFOREST (tReAtment benefIt prediction using raNdom FOREST)

![RainForest](/docs/img/RainForest.png)

# Usage instructions

***Load dataset***
```
snps = pd.read_csv("data/simulatedSNPdata_Survival.csv",sep=";") 
labels =  pd.read_csv("data/simulatedPatientinfo_Survival.csv",sep=";")
X_train, X_test, y_train, y_test = train_test_split(snps.iloc[:600,:200], labels.iloc[:600,:], test_size=0.33, random_state=42)
X_train.set_index("ID"), y_train.set_index("ID"), X_test.set_index("ID")
```

**Create a new RainForest instance**
```
from classifiers.ensembl import RainForest
RF = RainForest(max_depth=2, min_size=30, n_trees=20, n_features=2000)
```

**Fit a Trainingset**
```
RF.fit(X_train, y_train)
```

**Predict Treatment outcome**
```
X_test.apply(RF.ensembl_predict, axis=1)
```

**Serialize/deserialize a trained RF model**
```
import pickle
pickle.dump(RF.ttrees,open( "RF_trees.bin", "wb" ))
```
***Create new instance and load pickled trees***
```
RF_new = RainForest(max_depth=2, min_size=30, n_trees=20, n_features=2000)
RF_new.ttrees = pickle.load( open( "RF_trees.bin", "rb" ))
X_test.apply(RF_new.ensembl_predict, axis=1)
```



