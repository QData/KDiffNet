evaluate <- function(real,prediction){
testGraph = abs(prediction) > 0
realtruth = abs(real) > 0
tPa = sum(testGraph & realtruth) 
tNa = sum((testGraph == 0) & (realtruth == 0)) 
fPa = sum((testGraph == 1) & (realtruth == 0))  
fNa = sum((testGraph == 0) & (realtruth == 1))
precisiona=(tPa)/(tPa+fPa)
recalla=(tPa)/(tPa+fNa)
f1scorea=(2*precisiona*recalla)/(precisiona+recalla)
if(is.nan(f1scorea)){
f1scorea=0.0
}
f1scorea
}




