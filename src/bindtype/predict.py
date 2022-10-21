import os
import pandas as pd
import sklearn_json as skljson

from .descriptors import getDescriptors

def BindingTypeClassifier(data = None, data_path : str = None, output_path : str = None, smiles_column : str = 'SMILES'):
    """Predict binding type (allo- or orthosteric) of given molecule-protein pairs."""

    if data_path:
        data = pd.read_csv(data_path, sep='\t')
    elif data is None:
        raise ValueError('No data provided. Provide either data or data_path')
    
    # Get descriptors
    pairs = [ (smiles, target) for smiles, target in zip(data[smiles_column], data['accession']) ]
    x = getDescriptors(
        pairs, 
        scaler_from_file=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'models/proteinDesciptorScaler.json'),
    )

    # Load model
    print('Loading model...')
    model = skljson.from_json(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'models/BindingTypeClassifier_ClassAGPCR.json'))
    print('Predicting binding type...')
    preds = model.predict(x)
    probs = model.predict_proba(x)

    data['binding_type_prob'] = probs[:,1]
    data['binding_type_pred'] = [ 'orthosteric' if p == 0 else 'allosteric' for p in preds]
    
    # Save predictions
    if output_path:
        print(f'Saving predictions to {output_path}...')
        data.to_csv(output_path, sep='\t', index=False)

    return data