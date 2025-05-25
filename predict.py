import pandas as pd
import numpy as np
import joblib
from ProteinDescriptors import toBLOSUM62, toCTD
#加载训练好的模型
std_scaler=joblib.load("AOP.StandardScaler.joblib")

model=joblib.load("SVM.FusedAOPP.test01.CTD.BLOSUM62.OptF5F300.results_SVMmodelWith80Features.joblib")
#获取最佳特征组合
topFeats=pd.read_csv("AOPP.Test01.FusedCTD.BLOSUM62_Train.lgbm.sf.csv",header=0)

def getFeats(inFasta):
    blosum62_feats=toBLOSUM62(inFasta)
    blosum62_feats = blosum62_feats.drop(columns=['PID'])
    ctd_feats=toCTD(inFasta)
    ctd_feats = ctd_feats.drop(columns=['PID'])
    combined_feats = pd.concat([
        blosum62_feats.reset_index(drop=True),
        ctd_feats.reset_index(drop=True)
    ], axis=1)
    print(combined_feats.shape)
    print(combined_feats.head())
    return combined_feats
feats=getFeats("AOPP.test.fasta")
feats=feats[topFeats.columns]
feats_std=std_scaler.transform(feats)
pred_labels=model.predict(feats_std[:,:80])
pred_probas=model.predict_proba(feats_std[:,:80])
pred_results=pd.DataFrame()

pred_results["pred_label"]=pred_labels
pred_results["pred_proba"]=pred_probas[:,1]
pred_results.to_csv("测试_预测结果.csv")
