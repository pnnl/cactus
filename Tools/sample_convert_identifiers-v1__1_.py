#!/usr/bin/env python
# coding: utf-8

# In[3]:


get_ipython().system('pip install pubchempy')
get_ipython().system('pip install chembl_webresource_client')
get_ipython().system('pip install smilite')


# In[7]:


import pubchempy as pcp
from chembl_webresource_client.new_client import new_client
import smilite


# In[8]:


def common_name_to_SMILES(input_name: str)-> str:

    if(isinstance(input_name,str)):
        result=pcp.get_compounds(str(input_name),'name')[0]
        #print (result)
        return result.isomeric_smiles
    else:
        raise ValueError("Invalid input")
    


# In[9]:


def name_to_SMILES(input_name: str)-> str:
#checking if the input is a string variable or not
    if(isinstance(input_name,str)):
#checking for the first compund with the given input name
        result=pcp.get_compounds(str(input_name),'name')[0]
        #print (result)
        return result.isomeric_smiles
    else:
#if the input is not a string error will be raised
        raise ValueError("Invalid input")
    


# In[10]:


def cas_to_SMILES(input_query: str)-> str:

    if(isinstance(input_query,str)):
         result=pcp.get_compounds(str(input_query),'name')[0]
         #print (result)
         return result
#         get_cid("input_query", from = "xref/RN")
#         return pc_sect(702,"canonical smiles")
    else:
        raise ValueError("Invalid input")
    


# In[ ]:


def chemblid_to_smiles(input_id: str)-> str:
    
    if(isinstance(input_id,str)):
        from chembl_webresource_client.new_client import new_client
        molecule = new_client.molecule
        molecule.set_format('json')
        record_via_client = molecule.get(input_id)
        smiles_from_json = record_via_client['molecule_structures']['canonical_smiles']
        return smiles_from_json
    else:
        raise ValueError("Invalid input")
    


# In[11]:


def cid_to_SMILES(input_id: int)-> str:
#checking if the input is a integer variable or not
    if(isinstance(input_id,int)):
        result=pcp.Compound.from_cid(int(input_id))
        #print (result)
        return result.isomeric_smiles
    else:
#if the input is not a string error will be raised
        raise ValueError("Invalid input")
    


# In[12]:


def molecular_formula_to_SMILES(input_formula: str)-> str:
#checking if the input is a string variable or not
    if(isinstance(input_formula,str)):
#checking for the first compund with the given molecular formula
        result=pcp.get_compounds(str(input_formula),'formula')[0]
        #print (result)
        return result.isomeric_smiles
    else:
#if the input is not a string error will be raised
        raise ValueError("Invalid input")


# In[ ]:


def zinc_id_to_SMILES(input_id: str)-> str:
    #checking if the input is a string variable or not
    import smilite
    if(isinstance(input_id,str)):
#checking for the first compund with the given molecular formula
        smile_str = smilite.get_zinc_smile(input_id, backend='zinc15')
        #print (result)
        return smile_str
    else:
#if the input is not a string error will be raised
        raise ValueError("Invalid input")


# In[20]:


#i1=int(input("Enter the query "))

i2=input("Enter the query ")
#k=common_name_to_SMILES(i1)
#l=iupac_name_to_SMILES(i1)
m=common_name_to_SMILES(i2)
#o=cid_to_SMILES(i1)
#p=molecular_formula_to_SMILES(i2)
print(m)


# In[ ]:




