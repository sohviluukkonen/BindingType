import os
import pandas as pd
import sklearn_json as skljson

from collections import Counter
from imblearn.over_sampling import SMOTE 

from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.metrics import matthews_corrcoef, accuracy_score, precision_score, recall_score, confusion_matrix

from .descriptors import getDescriptors

def validate_model(data_path : str):

    # Load data
    data = pd.read_csv(data_path, sep='\t')
    data = data.drop_duplicates()

    # Get descriptors for training set
    print('Preparing training data...')
    train = data[data['Subset'] == 'train']
    train_pairs = [ (smiles, target) for smiles, target in zip(train['canonical_smiles'], train['accession']) ]
    train_x = getDescriptors(train_pairs, scaler_to_file='models/proteinDesciptorScaler.json')
    train_y = train['binding_binary'].values
    print(f'Original number of training data points per class: {Counter(train_y)}')
    sm = SMOTE(random_state=42)
    train_x, train_y = sm.fit_resample(train_x, train_y)
    print(f'Number of training data points per class after oversampling: {Counter(train_y)}')


    # Get descriptors for test set
    print('Preparing test data...')
    test = data[data['Subset'] == 'test']
    test_pairs = [ (smiles, target) for smiles, target in zip(test['canonical_smiles'], test['accession']) ]
    test_x = getDescriptors(test_pairs, scaler_from_file='models/proteinDesciptorScaler.json')
    test_y = test['binding_binary'].values

    # Train model
    print('Building and fitting the model...')
    model = RandomForestClassifier(n_estimators=1000, max_depth=10, random_state=42, n_jobs=64, verbose=True)
    model.fit(train_x, train_y)

    # Evaluate model
    print('Evaluating the model...')
    train_prob = model.predict_proba(train_x)[:,1]
    test_prob = model.predict_proba(test_x)[:,1]
    train_auc = roc_auc_score(train_y, train_prob)
    test_auc = roc_auc_score(test_y, test_prob)

    train_pred = model.predict(train_x)
    test_pred = model.predict(test_x)
    train_mcc = matthews_corrcoef(train_y, train_pred)
    test_mcc = matthews_corrcoef(test_y, test_pred)
    train_acc = accuracy_score(train_y, train_pred)
    test_acc = accuracy_score(test_y, test_pred)
    # train_pre = precision_score(train_y, train_pred)
    # test_pre = precision_score(test_y, test_pred)
    train_sen = recall_score(train_y, train_pred) 
    test_sen = recall_score(test_y, test_pred)
    tn, fp, fn, tp = confusion_matrix(train_y, model.predict(train_x)).ravel()
    train_spe = tn / (tn+fp)
    tn, fp, fn, tp = confusion_matrix(test_y, model.predict(test_x)).ravel()
    test_spe = tn / (tn+fp)

    print('----------------------------------------')
    print('Train AUC: {:.3f}'.format(train_auc))
    print('Train MCC: {:.3f}'.format(train_mcc))
    print('Train ACC: {:.3f}'.format(train_acc))
    #print('Train PRE: {:.3f}'.format(train_pre))
    print('Train SEN: {:.3f}'.format(train_sen))
    print('Train SPE: {:.3f}'.format(train_spe))
    print('----------------------------------------')
    print('Test AUC: {:.3f}'.format(test_auc))
    print('Test MCC: {:.3f}'.format(test_mcc))
    print('Test ACC: {:.3f}'.format(test_acc))
    #print('Test PRE: {:.3f}'.format(test_pre))
    print('Test SEN: {:.3f}'.format(test_sen))
    print('Test SPE: {:.3f}'.format(test_spe))
    print('----------------------------------------')

def train_model(data_path : str, model_name : str):

    # Get descriptors for training set
    print('Preparing training data...')
    train =  pd.read_csv(data_path, sep='\t')
    train_pairs = [ (smiles, target) for smiles, target in zip(train['canonical_smiles'], train['accession']) ]
    train_x = getDescriptors(train_pairs, scaler_to_file='models/proteinDesciptorScaler.json')
    train_y = train['binding_binary'].values
    print(f'Original number of training data points per class: {Counter(train_y)}')
    sm = SMOTE(random_state=42)
    train_x, train_y = sm.fit_resample(train_x, train_y)
    print(f'Number of training data points per class after oversampling: {Counter(train_y)}')

    # Train model
    print('Building and fitting the model...')
    model = RandomForestClassifier(n_estimators=1000, max_depth=10, random_state=42, n_jobs=64, verbose=True)
    model.fit(train_x, train_y)

    # Save model
    print('Saving the model...')
    skljson.to_json(model, f'models/{model_name}')

if __name__ == '__main__':

    os.makedirs('models', exist_ok=True)

    # Build and evaluate model binding type classification and train prodution model
    validate_model('../../data/BindingTypeDataset_ClassAGPCRs.tsv')
    train_model('../../data/BindingTypeDataset_ClassAGPCRs.tsv', 'BindingTypeClassifier_ClassAGPCR.json')
