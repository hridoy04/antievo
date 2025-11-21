# %%
#Load the libraries

#!/usr/bin/env python3
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import numpy as np
import torch
import esm
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from numpy import dot
from numpy.linalg import norm
from Bio import SeqIO
from scipy.special import softmax
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm
import bz2
import pickle
import _pickle as cPickle
import re
import matplotlib.pyplot as plt
import argparse,csv,sys
import os,requests
import torch.nn as nn
import torch.nn.functional as nnF
from models import *
import glob
from typing import List, Dict


# %%
model_path=lambda x:'models/model_'+str(x)+'.pts'
github_url=lambda x:f"https://github.com/ComputBiophys/ProtRAP-LM/releases/download/Version1.0/model_{str(x)}.pts"

def download_file(url, output_path):
    try:
        response = requests.get(url)
        response.raise_for_status()  # Check whether the request was successful
        with open(output_path, 'wb') as f:
            f.write(response.content)
        print(f"Downloaded file from {url} to {output_path}")
    except Exception as e:
        print(f"Error downloading file: {e}, You may manually download this one")

for i in range(10):
    if not os.path.exists(model_path(i)):
        print('Downloading model_'+str(i))
        download_file(github_url(i), model_path(i))


def fasta_load(fasta_dir):
    fp = open(fasta_dir, 'r')
    lines = fp.readlines()
    fp.close()
    sequence = ''
    for line in lines[1:]:
        sequence = sequence + line.split()[0]
    return sequence
def weight_MSE_loss(labels,logits,weights=1):
    l=(labels-logits)**2
    l=l*weights
    return torch.sum(l)
def focal_loss_softmax(labels,logits):
    y_pred=logits
    l=-labels*torch.log(y_pred+1e-8)*((1-y_pred)**2)
    return torch.sum(l)

class MultiScaleCNN(nn.Module):
    def __init__(self,input_dim=1280,output_dim=256):#,size=[3,7,11],padding=[1,3,5]):
        super().__init__()
        self.cnn1=nn.Conv1d(input_dim,output_dim,3,padding=1)
        self.cnn2=nn.Conv1d(input_dim,output_dim,5,padding=2)
        self.cnn3=nn.Conv1d(input_dim,output_dim,7,padding=3)
        self.cnn4=nn.Conv1d(input_dim,output_dim,9,padding=4)
    def forward(self,x):
        x=x.permute(0,2,1)
        x1=self.cnn1(x)
        x2=self.cnn2(x)
        x3=self.cnn3(x)
        x4=self.cnn4(x)
        x=torch.cat((x1,x2,x3,x4), -2)
        x=x.permute(0,2,1)
        return x
        
class ProtRAP_LM_Model(nn.Module):
    def __init__(self,input_dim=1280,n_hidden=256,num_layers=2,dropout=0.1):
        super().__init__()
        assert n_hidden%8==0

        self.keep_prob=1-dropout
        self.begin_linears=nn.Sequential(
            nn.Linear(input_dim,n_hidden*2),nn.ReLU(),nn.Dropout(self.keep_prob),
            nn.Linear(n_hidden*2,n_hidden),)
        self.cnn=MultiScaleCNN(input_dim=n_hidden,output_dim=int(n_hidden/4))
        encoder_layer=nn.TransformerEncoderLayer(d_model=n_hidden, nhead=4,activation='gelu',batch_first=True)
        self.encoder= nn.TransformerEncoder(encoder_layer,num_layers=num_layers)
        self.pred=nn.Sequential(
            nn.Linear(n_hidden,int(n_hidden/2)),nn.ReLU(),nn.Dropout(self.keep_prob),nn.Linear(int(n_hidden/2),64),
            nn.Linear(64,3),nn.Sigmoid())
        return
    def forward(self,x):
        x=self.begin_linears(x)
        x=self.cnn(x)+x
        x=self.encoder(x)
        prediction=self.pred(x)
        return prediction

class ProtRAP_LM():

    def __init__(self,device_name='cpu'):
        device = torch.device(device_name)
        self.device=device

        esm_model, alphabet = torch.hub.load("facebookresearch/esm:main", "esm2_t33_650M_UR50D")
        batch_converter = alphabet.get_batch_converter()
        esm_model=esm_model.eval().to(device)
        models=[]
        for i in range(10):
            model=torch.jit.load('./models/model_'+str(i)+'.pts').to(device).eval()
            models.append(model)
        self.models=models
        self.esm_model=esm_model
        self.batch_converter=batch_converter
        
    def predict(self,seq):
        data=[('prot',seq)]
        _, _, batch_tokens = self.batch_converter(data)
        batch_tokens=batch_tokens.to(self.device)
        preds=[]
        with torch.no_grad():
            results=self.esm_model(batch_tokens,repr_layers=[33])
            Repr= results["representations"][33]
            for model in self.models:
                pred=model(Repr).to(torch.device("cpu"))
                preds.append(np.array(pred[0,1:-1,:]))
        preds=np.array(preds)
        mean_pred=np.mean(preds,axis=0)
        std_pred=np.std(preds,axis=0)
        result=np.concatenate((mean_pred,std_pred),axis=-1)
        return result
def run_ProtRAP_LM(fasta_file, result_rsa_dir, device):

    seq = fasta_load(fasta_file)
    result=ProtRAP_LM(device).predict(seq)
    np.savetxt(result_rsa_dir + 'results.txt', np.column_stack((result[:,0:1], result[:,1:2], result[:,0:1]*result[:,1:2], (1-result[:,0:1])*result[:,1:2], 1-result[:,1:2])), fmt='%.4f %.4f %.4f %.4f %.4f', header='MCP RASA RLA RSA RBSA')


# %%
#this part will run the model and save the RSA results in csv (plus txt) format and in dataframe that can be used for further analysis

fasta_file = input("Enter the path to the fasta file: ")
result_rsa_dir = input('name the directory to save the results:')

os.makedirs(result_rsa_dir, exist_ok=True)
device = 'cpu' #can use cuda:0 in gpu while running in gpu node or PC
if __name__=='__main__':
    
    run_ProtRAP_LM(str(fasta_file), str(result_rsa_dir), str(device))

rsa = pd.DataFrame({'SurfaceAccessibility' : np.loadtxt(result_rsa_dir + 'results.txt', comments="#", usecols=3)})
seq = [i.seq for i in SeqIO.parse(fasta_file, 'fasta')]
seq = str(seq[0])
id = [i.id for i in SeqIO.parse(fasta_file, 'fasta')]
id = str(id[0])
df_seq = pd.DataFrame({'position': list(range(1, len(seq)+1)), 'ref': list(seq)})
df_rsa = pd.concat([df_seq, rsa], axis=1)
print(df_rsa)
df_rsa.to_csv(result_rsa_dir + '/' + id + 'RSA.csv', index=False)


#alternate methods:
# pd.read_csv("rsaresults/results.txt", sep=r"\s+", comment="#", names=["MCP","RASA","RLA","RSA","RBSA"]) 
# df[["RSA"]]


"""
SEMA-1D inference - simple function for single sequence prediction.

Example usage:
    from sema1d_lib import predict_sequence
    
    # Your sequence variable
    sequence = "MKTIIALSYIFCLVFA"
    
    # Get predictions
    scores = predict_sequence(
        sequence=sequence,
        model_type="esm2_650m",
        weights_glob="../models/sema_1d_ESM2_*.pth"
    )
    
    # scores is a list of floats, one per amino acid
    print(scores)  # [0.52, 0.48, 0.61, ...]
"""


def predict_sequence(
    sequence: str,
    model_type: str,
    weights_glob: str,
    device: str = None,
    max_len: int = None
) -> List[float]:
    """
    Predict SEMA-1D scores for a single sequence.
    
    Args:
        sequence: Amino acid sequence string
        model_type: One of "esm1v", "esm2_650m", "esm2_3b"
        weights_glob: Glob pattern for checkpoints, e.g., "../models/*.pth"
        device: "cuda" or "cpu" (auto-detected if None)
        max_len: Maximum sequence length (e.g., 1022), None for no limit
        
    Returns:
        List of scores (floats), one per amino acid in the sequence
    """
    device = torch.device(device if device else ("cuda" if torch.cuda.is_available() else "cpu"))
    
    # Truncate if needed
    if max_len is not None:
        sequence = sequence[:max_len]
    
    # Build model
    model, alphabet, rep_layer = _build_model(model_type)
    batch_converter = alphabet.get_batch_converter()
    
    # Load checkpoints
    ckpts = sorted(glob.glob(weights_glob))
    if not ckpts:
        raise ValueError(f"No checkpoints matched: {weights_glob}")
    
    ensembles = []
    for ck in ckpts:
        m, _, _ = _build_model(model_type)
        m.to(device).eval()
        _load_weights(m, ck)
        ensembles.append(m)
    
    print(f"[INFO] Loaded {len(ensembles)} checkpoints")
    
    # Tokenize
    _, _, toks = batch_converter([("seq", sequence)])
    toks = toks.to(device)
    
    # Predict with ensemble
    with torch.no_grad():
        preds = []
        for m in ensembles:
            rep = m["backbone"].forward(toks, repr_layers=[rep_layer])["representations"][rep_layer]
            rep = rep[:, 1:-1, :]  # drop BOS/EOS
            logits = m["head"](rep)[:, :, 1]  # class-1 logit
            preds.append(logits.squeeze(0).float().cpu())
        
        # Average predictions
        scores = torch.stack(preds, 0).mean(0).tolist()
    
    return scores


def _build_model(kind: str):
    """Build model architecture."""
    if kind == "esm1v":
        backbone, alphabet = esm.pretrained.esm1v_t33_650M_UR90S_1()
        rep_layer = 33
    elif kind == "esm2_650m":
        backbone, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
        rep_layer = 33
    elif kind == "esm2_3b":
        backbone, alphabet = esm.pretrained.esm2_t36_3B_UR50D()
        rep_layer = 36
    else:
        raise ValueError("kind must be: esm1v, esm2_650m, or esm2_3b")
    
    embed_dim = getattr(backbone, "embed_dim", 1280)
    head = nn.Linear(embed_dim, 2)
    model = nn.ModuleDict({"backbone": backbone, "head": head})
    return model, alphabet, rep_layer


def _load_weights(model: nn.Module, ckpt_path: str):
    """Load checkpoint into model."""
    sd = torch.load(ckpt_path, map_location="cpu")
    if isinstance(sd, dict) and "state_dict" in sd:
        sd = sd["state_dict"]
    
    # Normalize keys
    mapped = {}
    for k, v in sd.items():
        if k.startswith("esm1v."):
            mapped["backbone." + k[6:]] = v
        elif k.startswith("classifier."):
            mapped["head." + k[11:]] = v
        else:
            mapped[k] = v
    
    model.load_state_dict(mapped, strict=False)


# Example usage
if __name__ == "__main__":
    # Your sequence from wherever (variable, file, database, etc.)
    my_sequence = seq
    
    # Get predictions
    scores = predict_sequence(
        sequence=my_sequence,
        model_type="esm2_3b",
        weights_glob="sema_1d_ESM2_0_3B.pth"
    )
    
    # Print results
    print(f"Sequence: {my_sequence}")
    print(f"Scores: {scores}")
    print(f"\nPer-position:")
    for i, (aa, score) in enumerate(zip(my_sequence, scores), 1):
        print(f"  Position {i}: {aa} -> {score:.3f}")

#saving the antigenicity results in dataframe and csv

antigenicity = pd.DataFrame({'Antigenicity' : scores})
antigenicity.to_csv(result_rsa_dir + '/' + id + 'Antigenicity.csv', index=False)
df_rsa_antigenicity = pd.concat([df_rsa, antigenicity], axis=1)
print(df_rsa_antigenicity)



#functions to decompress pickle files and zero matching index columns
def decompress_pickle(file):
  data = bz2.BZ2File(file, 'rb')
  data = cPickle.load(data)
  return data

def zero_matching_index_column(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    
    for idx in df.index:
        match = re.search(r'(\d+)', idx)
        if match:
            col = int(match.group(1))
            if col in df.columns:
                df.at[idx, col] = 0  # Only set the matching column to 0
    return df



# %%
############################################################################################################
######################################### Imports ##########################################################
######################################### Pytorch imports ################################################## 
import torch
from torch.utils.data import DataLoader
######################################### Transformers imports ############################################## 
# Use a pipeline as a high-level helper
from transformers import pipeline
# Load model directly
from transformers import AutoTokenizer, AutoModelForMaskedLM
######################################### Numpy and Pandas imports ############################################
import numpy as np
import pandas as pd
######################################### Biopython imports ##################################################
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Entrez
from Bio import SeqIO
######################################### SciPy imports ##################################################
#from scipy.special import softmax
######################################### tqdm imports ##################################################
from tqdm import tqdm
############################################################################################################
####################################### Functions ##########################################################
from HuggingFace_Functions import *
import os
############################################################################################################
################################ Input parameters ##########################################################
params = {
    #Specify ORF to analyse
    'specify_orf':False,
    #Decide if <mask> is used in DMS
    'include_mask_token_in_DMS' : True,
    #Select model
    'model_name' : "facebook/esm2_t33_650M_UR50D",
    #Specify if ids are from genbank or are protein names
    'use_genbank':False,
    #Specify top directory to save DMS
    'container_directory':id,
    #Specify if genbank is a virus assembly (currently just influenza viruses)
    'ncbi_virus_assembly':False,
    'save_representations':True
}


############################################################################################################
################################ Global Variables ##########################################################
#Amino acid tokens
amino_acids = ["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
############################################################################################################
################################ Load Model into GPU #######################################################
#Load model and tokenizer
model = AutoModelForMaskedLM.from_pretrained(params['model_name'],output_hidden_states = True)
tokenizer = AutoTokenizer.from_pretrained(params['model_name'])

#Assign device for using model (GPU)
device = torch.device("cuda:0")
model = model.to(device)
################################################################################################################
############################################## Main ############################################################

#


id_list = {
    id: seq
}

if params['ncbi_virus_assembly'] == True:
    assembly = []

for identifier,query in id_list.items():
    #Create directory for protein/genbank file
    os.makedirs(params['container_directory']+'/'+identifier, exist_ok=True)
    if params['save_representations'] == True:
        os.makedirs(params['container_directory']+'/'+identifier+'/representations', exist_ok=True)
    #Use genbank file for DMS or just single sequence
    if params['use_genbank'] == True:
        #Retrieve genbank file
        Entrez.email = "sample@example.org"
        handle = Entrez.efetch(db="nucleotide",id=query,rettype="gb",retmode="gb")
        genbank_file = SeqIO.read(handle, "genbank")

        #Translate nucleotide to proteins using genbank
        Coding_Regions= translate_with_genbank(genbank_file.seq,genbank_file)
        print(Coding_Regions)
        Mature_Proteins= translate_mat_proteins_with_genbank(genbank_file.seq,genbank_file)
        if len(Mature_Proteins) != 0:
            polyprotein_orfs =set([Coding_Regions[prot]["ORF"] for prot in Coding_Regions.keys()])
            polyprotein_orfs =set([Mature_Proteins[prot]["ORF"] for prot in Mature_Proteins.keys()])
            Filtered_Coding_Regions = {**Coding_Regions}
            for orf in Coding_Regions.keys():
                if Coding_Regions[orf]["ORF"] in polyprotein_orfs:
                    del Filtered_Coding_Regions[orf]
            Merged_Coding_Regions = {**Filtered_Coding_Regions,**Mature_Proteins}
        else:
            Merged_Coding_Regions = Coding_Regions
        embeddings = {}
        if params['specify_orf'] !="":
            Merged_Coding_Regions = {params['specify_orf']:Merged_Coding_Regions[params['specify_orf']]}
    else:
        Merged_Coding_Regions = {identifier:{"Sequence":query}}

    ############################################################################################################
    ############################################## Inner Main ########################################################
    all_dfs = []

    

    
    #Embed Sequence
    for key,value in Merged_Coding_Regions.items():

        #Create dictionary for logits and embeddings
        representations = {}
        representations[key] = {}
        print(value['Sequence'])
        
        dfs = []
        # Define the reference sequence
        reference_sequence = value['Sequence']

        #Define sequence length
        sequence_length = len(reference_sequence)

        #Perform and batch DMS on reference sequence
        ids,sequences = Batch_DMS(reference_sequence,mask=params['include_mask_token_in_DMS'])

        #Calculate reference logits and embedding
        reference_logits,reference_embedding,reference_tokens = embed_batch([reference_sequence],tokenizer,model,device)
        reference_logits,reference_embedding,reference_tokens =reference_logits[0],reference_embedding[0],reference_tokens[0]
        reference_sequence_grammaticality = np.sum(reference_logits[np.arange(reference_logits.shape[0]), reference_tokens])


        #Embed sequences in batches (one batch is equivelent to all the amino acids at a position in the sequence)
        for batch_number in tqdm(range(0,len(sequences))):
            #Index batch
            sequence_zero_position = batch_number

            batch_ids = ids[batch_number]
            batch_sequences = sequences[batch_number]

            #Embed batch and extract logits and mean embeddings
            logits,embeddings,tokens = embed_batch(batch_sequences,tokenizer,model,device)

            #Calculate semantic scores
            semantic_score = abs(embeddings - reference_embedding).sum(axis=1)
            
            logit_index =logits[:,sequence_zero_position,:]

            #Retrieve reference grammaticality (Log-likelihood of reference amino acid)
            reference_grammaticality = reference_logits[sequence_zero_position,reference_tokens[sequence_zero_position]]

            #Retrieve grammaticality (Log-likelihood of mutant amino acid)
            grammaticality = reference_logits[sequence_zero_position,tokens[:,sequence_zero_position]]
           
            #Calculate relative grammaticality (LLR)
            relative_grammaticality = grammaticality - reference_grammaticality

            sequence_logits = np.array([logits[:,i,:][np.arange(logit_index.shape[0]), tokens[:,i]]  for i in range(0,sequence_length)])

            #Calculate sequence grammaticality
            sequence_grammaticality = np.sum(sequence_logits,axis=0)

            #Calculate relative sequence grammaticality (PLLR)
            relative_sequence_grammaticality = sequence_grammaticality - reference_sequence_grammaticality

            #Calculate mutated_grammaticality
            mutated_grammaticality = logit_index[np.arange(logit_index.shape[0]), tokens[:,sequence_zero_position]]

            #Calculate relative_mutated_grammaticality (Reference sequence logits for sequence minus mutated sequence logits for sequence)
            relative_mutated_grammaticality = mutated_grammaticality - reference_grammaticality
            
            #<mask> grammaticalities
            if params['include_mask_token_in_DMS'] == True:
                #<mask> token is the last sequence in the DMS batch
                masked_grammaticality = logit_index[logit_index.shape[0]-1,tokens[:,sequence_zero_position]]
                relative_masked_grammaticality = masked_grammaticality - reference_grammaticality
            else:
                masked_grammaticality = np.arange(20,np.nan)
                relative_masked_grammaticality = np.arange(20,np.nan)

            print(sequence_logits.shape)
            #Save Representations (embeddings and Logits)
            if params['save_representations'] == True:
                for i in range(0,len(batch_ids)):
                    representations[key][batch_ids[i]] = {}
                    representations[key][batch_ids[i]]["logits"] = sequence_logits[:,i]
                    representations[key][batch_ids[i]]["embeddings"] = embeddings[i]
    
            #Append to dataframe
            df = pd.DataFrame([batch_ids, semantic_score, relative_grammaticality, relative_sequence_grammaticality,relative_mutated_grammaticality,relative_masked_grammaticality,grammaticality,sequence_grammaticality,mutated_grammaticality,masked_grammaticality,np.exp(grammaticality),np.exp(mutated_grammaticality),np.exp(masked_grammaticality)],
                            index = ['label','semantic_score','relative_grammaticality','relative_sequence_grammaticality','relative_mutated_grammaticality','relative_masked_grammaticality','grammaticality','sequence_grammaticality','mutated_grammaticality','masked_grammaticality','probability','mutated_probability','masked_probability']).T
            df['ref'] = df.label.str[0]
            df['alt'] = ['<mask>' if '<mask>' in i else i[-1] for i in  df.label]
            df['position'] = [int(i[1:-len('<mask>')]) if '<mask>' in i else int(i[1:-1]) for i in  df.label]
            df['region'] = key
            df['subunit'] = ""
            df['domain']  = ""
            df['GenbankID_or_Pro']  = identifier
            df['reference_grammaticality'] = reference_grammaticality
            dfs.append(df)
        dms_table = pd.concat(dfs,axis=0)
        dms_table = dms_table.sort_values('semantic_score')
        dms_table['semantic_rank'] = dms_table.reset_index().index.astype(int) + 1
        dms_table = dms_table.sort_values('grammaticality')
        dms_table['grammatical_rank'] =dms_table .reset_index().index.astype(int) + 1
        dms_table['acquisition_priority'] = dms_table['semantic_rank'] + dms_table['grammatical_rank']

        dms_table = dms_table.sort_values('sequence_grammaticality')
        dms_table['sequence_grammatical_rank'] =dms_table.reset_index().index.astype(int) + 1
        dms_table['sequence_acquisition_priority'] = dms_table['semantic_rank'] + dms_table['sequence_grammatical_rank']
        dms_table.to_csv(params['container_directory']+'/'+identifier+'/'+'DMS_'+key+'.csv',index=False)
        all_dfs.append(dms_table)
        compressed_pickle(params['container_directory']+'/'+identifier+'/representations/DMS_'+key+'.pbz',representations)


    ############################################################################################################
    ################################ Output DMS as DataFrame ###################################################
    all_dfs =  pd.concat(all_dfs ,axis=0)
    all_dfs.to_csv(params['container_directory']+'/'+identifier+"_DMS_scores.csv",index=False)
    if params['ncbi_virus_assembly'] == True:
        assembly.append(all_dfs)
if params['ncbi_virus_assembly'] == True:
    assembly = pd.concat(assembly,axis=0)
    assembly.to_csv(params['container_directory']+'/'+params['container_directory']+"_Assembly_DMS_scores.csv",index=False)


# %%
pbz2file = os.getcwd() + '/' + params['container_directory']+'/'+identifier+'/representations/DMS_'+key+'.pbz.pbz2'
logits_and_embeddings = decompress_pickle(pbz2file)
#print(logits_and_embeddings[id].keys())
first_match = [x for x in [i for i in logits_and_embeddings[id].keys()][0:20] if x[-1] == x[0]][0]
#print(first_match)
index = [k for k in logits_and_embeddings[id].keys()]
values_df = pd.DataFrame([logits_and_embeddings[id][v]['logits'] for v in logits_and_embeddings[id].keys()],index= index)
values_df.columns = [int(c)+1 if c !='label' else 'label' for c in values_df.columns]
values_df.index.name = 'label'
values_df = values_df.drop([i for i in values_df.index if '<mask>' in i])
values_df = np.exp(values_df) - np.exp(logits_and_embeddings[id][first_match]['logits'])
values_df = zero_matching_index_column(values_df)
values_df['position'] = [int(i[1:-1]) for i in values_df.index]
values_df = values_df[(values_df.position.isin(df_rsa_antigenicity.position))].drop(columns=['position'])
mean = np.ravel(values_df.values).mean()
std =  np.ravel(values_df.values).std()
significant_positions = (values_df[(values_df>mean+(std)) | (values_df<mean-(std)) ].isna() == False).apply(np.flatnonzero, axis=1) 
very_significant_positions = (values_df[(values_df>mean+(std*2)) | (values_df<mean-(std*2)) ].isna() == False).apply(np.flatnonzero, axis=1)

# %%
ranked_accessibilities = []
for i in range(len(significant_positions)):
    
    label = significant_positions.index[i]

    # Check if either significant_positions OR very_significant_positions has empty list
    if len(significant_positions.loc[label]) == 0 or len(very_significant_positions.loc[label]) == 0:
        # Use the mutation's own position as fallback for any empty categories
        mutation_position = int(significant_positions.index[i][1:-1])
        
        # Get accessibility and antigenicity values for the mutation's own position
        fallback_rsa = df_rsa_antigenicity.iloc[mutation_position - 1].SurfaceAccessibility
        fallback_antigenicity = df_rsa_antigenicity.iloc[mutation_position - 1].Antigenicity

        ranked_accessibilities.append(
            pd.Series({
             'label': significant_positions.index[i],
             'position': int(significant_positions.index[i][1:-1]),
             'ref': str(significant_positions.index[i][0]),
             'alt': str(significant_positions.index[i][-1]),

             # For any empty positions, use the fallback value directly
             'significant_accessibility_median': fallback_rsa,
             'very_significant_accessibility_median': fallback_rsa,
             'significant_accessibility_mean': fallback_rsa,
             'very_significant_accessibility_mean': fallback_rsa,
             'significant_accessibility_sum': fallback_rsa,
             'very_significant_accessibility_sum': fallback_rsa,
             'significant_accessibility_max': fallback_rsa,
             'very_significant_accessibility_max': fallback_rsa,
             
             'significant_antigenicity_median': fallback_antigenicity,
             'very_significant_antigenicity_median': fallback_antigenicity,
             'significant_antigenicity_mean': fallback_antigenicity,
             'very_significant_antigenicity_mean': fallback_antigenicity,
             'significant_antigenicity_sum': fallback_antigenicity,
             'very_significant_antigenicity_sum': fallback_antigenicity,
             'significant_antigenicity_max': fallback_antigenicity,
             'very_significant_antigenicity_max': fallback_antigenicity
            })
        )
    else:
        # Normal processing for mutations with significant positions
        significant_likelihoods_rsa = df_rsa_antigenicity[df_rsa_antigenicity.position.isin(significant_positions.loc[label])].SurfaceAccessibility
        very_significant_likelihoods_rsa = df_rsa_antigenicity[df_rsa_antigenicity.position.isin(very_significant_positions.loc[label])].SurfaceAccessibility 
        

        significant_likelihoods_antigenicity = df_rsa_antigenicity[df_rsa_antigenicity.position.isin(significant_positions.loc[label])].Antigenicity
        very_significant_likelihoods_antigenicity = df_rsa_antigenicity[df_rsa_antigenicity.position.isin(very_significant_positions.loc[label])].Antigenicity

        ranked_accessibilities.append(
            pd.Series({
             'label': significant_positions.index[i],
             'position': int(significant_positions.index[i][1:-1]),
             'ref': str(significant_positions.index[i][0]),
             'alt': str(significant_positions.index[i][-1]),

             # RSA metrics - no error handling needed since we filtered empty cases above
             'significant_accessibility_median': np.nanmedian(np.ravel(significant_likelihoods_rsa.to_numpy(dtype="float64", na_value=np.nan))),
             'very_significant_accessibility_median': np.nanmedian(np.ravel(very_significant_likelihoods_rsa.to_numpy(dtype="float64", na_value=np.nan))),
             'significant_accessibility_mean': np.nanmean(np.ravel(significant_likelihoods_rsa.to_numpy(dtype="float64", na_value=np.nan))),
             'very_significant_accessibility_mean': np.nanmean(np.ravel(very_significant_likelihoods_rsa.to_numpy(dtype="float64", na_value=np.nan))),
             'significant_accessibility_sum': np.nansum(np.ravel(significant_likelihoods_rsa.to_numpy(dtype="float64", na_value=np.nan))),
             'very_significant_accessibility_sum': np.nansum(np.ravel(very_significant_likelihoods_rsa.to_numpy(dtype="float64", na_value=np.nan))),
             'significant_accessibility_max': np.nanmax(np.ravel(significant_likelihoods_rsa.to_numpy(dtype="float64", na_value=np.nan))),
             'very_significant_accessibility_max': np.nanmax(np.ravel(very_significant_likelihoods_rsa.to_numpy(dtype="float64", na_value=np.nan))),
             
             # Antigenicity metrics - no error handling needed since we filtered empty cases above
             'significant_antigenicity_median': np.nanmedian(np.ravel(significant_likelihoods_antigenicity.to_numpy(dtype="float64", na_value=np.nan))),
             'very_significant_antigenicity_median': np.nanmedian(np.ravel(very_significant_likelihoods_antigenicity.to_numpy(dtype="float64", na_value=np.nan))),
             'significant_antigenicity_mean': np.nanmean(np.ravel(significant_likelihoods_antigenicity.to_numpy(dtype="float64", na_value=np.nan))),
             'very_significant_antigenicity_mean': np.nanmean(np.ravel(very_significant_likelihoods_antigenicity.to_numpy(dtype="float64", na_value=np.nan))),
             'significant_antigenicity_sum': np.nansum(np.ravel(significant_likelihoods_antigenicity.to_numpy(dtype="float64", na_value=np.nan))),
             'very_significant_antigenicity_sum': np.nansum(np.ravel(very_significant_likelihoods_antigenicity.to_numpy(dtype="float64", na_value=np.nan))),
             'significant_antigenicity_max': np.nanmax(np.ravel(significant_likelihoods_antigenicity.to_numpy(dtype="float64", na_value=np.nan))),
             'very_significant_antigenicity_max': np.nanmax(np.ravel(very_significant_likelihoods_antigenicity.to_numpy(dtype="float64", na_value=np.nan)))
            })
        )

# %%
ranked_accessibilities = pd.DataFrame(ranked_accessibilities)
ranked_accessibilities['virus_name'] = id.split('_')[0]
ranked_accessibilities['protein_name'] = id.split('_')[1]
ranked_accessibilities.to_csv(id+'_epistatic_accessibility_antigenicity_v2.csv')
ranked_accessibilities

print(ranked_accessibilities.head(n=50))

# %%
# create and save lineplot without displaying
'''fig, ax = plt.subplots(figsize=(12, 6))
ranked_accessibilities.plot(x='position', y='significant_accessibility_median', ax=ax)
ax.set_title(f'Epistatic surface accessibility of every possible mutation from the {id.split("_")[0]} {id.split("_")[1]}')
ax.set_xlabel('Position')
ax.set_ylabel('Significant Accessibility Median')
fig.savefig(f'{id}_epiacc_lineplot.png', dpi=300, bbox_inches='tight')
plt.close(fig)

# create and save scatterplot without displaying
fig, ax = plt.subplots(figsize=(20, 10))
ranked_accessibilities.plot(x='position', y='significant_accessibility_median', kind='scatter', ax=ax)
ax.set_title(f'Epistatic surface accessibility of every possible mutation from the {id.split("_")[0]} {id.split("_")[1]}')
ax.set_xlabel('Position')
ax.set_ylabel('Significant Accessibility Median')
fig.savefig(f'{id}_epiacc_scatterplot.png', dpi=300, bbox_inches='tight')
plt.close(fig)


# %%
# Create pivot table for heatmap
pivot_data = ranked_accessibilities.pivot_table(
    values='significant_accessibility_median',
    index='alt',
    columns='position',
    fill_value=0
)

fig, ax = plt.subplots(figsize=(15, 8))
sns.heatmap(pivot_data,
            annot=False,
            cmap='Reds',
            cbar_kws={'label': 'Significant Accessibility Median'},
            ax=ax)
ax.set_title('Epistatic Accessibility by Position and Amino Acid')
ax.set_xlabel('Position')
ax.set_ylabel('Alternative Amino Acid')
fig.tight_layout()
fig.savefig(f'{id}_epiacc_heatmap.png', dpi=300, bbox_inches='tight')
plt.close(fig)

# create and save lineplot for antigenicity without displaying
fig, ax = plt.subplots(figsize=(12, 6))
ranked_accessibilities.plot(x='position', y='significant_antigenicity_median', ax=ax)
ax.set_title(f'Epistatic antigenicity of every possible mutation from the {id.split("_")[0]} {id.split("_")[1]}')
ax.set_xlabel('Position')
ax.set_ylabel('Significant Antigenicity Median')
fig.savefig(f'{id}_epiant_lineplot.png', dpi=300, bbox_inches='tight')
plt.close(fig)

# create and save scatterplot for antigenicity without displaying
fig, ax = plt.subplots(figsize=(20, 10))
ranked_accessibilities.plot(x='position', y='significant_antigenicity_median', kind='scatter', ax=ax)
ax.set_title(f'Epistatic antigenicity of every possible mutation from the {id.split("_")[0]} {id.split("_")[1]}')
ax.set_xlabel('Position')
ax.set_ylabel('Significant Antigenicity Median')
fig.savefig(f'{id}_epiant_scatterplot.png', dpi=300, bbox_inches='tight')
plt.close(fig)

# create and save heatmap for antigenicity
pivot_data_ant = ranked_accessibilities.pivot_table(
    values='significant_antigenicity_median',
    index='alt',
    columns='position',
    fill_value=0
)

fig, ax = plt.subplots(figsize=(15, 8))
sns.heatmap(pivot_data_ant,
            annot=False,
            cmap='Reds',
            cbar_kws={'label': 'Significant Antigenicity Median'},
            ax=ax)
ax.set_title('Epistatic Antigenicity by Position and Amino Acid')
ax.set_xlabel('Position')
ax.set_ylabel('Alternative Amino Acid')
fig.tight_layout()
fig.savefig(f'{id}_epiant_heatmap.png', dpi=300, bbox_inches='tight')
plt.close(fig)'''