#!/usr/bin/env python
# coding: utf-8

# In[3]:


get_ipython().system('pip install pubchempy')
get_ipython().system('pip install chembl_webresource_client')
get_ipython().system('pip install smilite')


# In[8]:


import pubchempy as pcp
from chembl_webresource_client.new_client import new_client
import smilite


# In[13]:


def inchikey_to_SMILES(input_name: str)-> str:
    """
    Convert inchikey into SMILES notation.

    Parameters:
    input_name (str): The InchIkey of the chemical compound.

    Returns:
    str: The SMILES notation in the output format.
    """
    if(isinstance(input_name,str)):
        result=pcp.get_compounds(str(input_name),'inchikey')[0]
        #print (result)
        return result.isomeric_smiles
    else:
        raise ValueError("Invalid input")
    


# In[16]:





# In[9]:


def name_to_SMILES(input_name: str)-> str:
    """
    Convert chemical name into SMILES notation.

    Parameters:
    input_name (str): The name of the chemical compound.

    Returns:
    str: The SMILES notation in the output format.
    """

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
    """
    Convert cas number into SMILES notation.

    Parameters:
    input_query (str): The cas number of the chemical compound.

    Returns:
    str: The SMILES notation in the output format.
    """
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
    """
    Convert ChEMBL id into SMILES notation.

    Parameters:
    input_id (str): The ChEMBL database id of the chemical compound.

    Returns:
    str: The SMILES notation in the output format.
    """    
    if(isinstance(input_id,str)):
        from chembl_webresource_client.new_client import new_client
        molecule = new_client.molecule
        molecule.set_format('json')
        record_via_client = molecule.get(input_id)
        smiles_from_json = record_via_client['molecule_structures']['canonical_smiles']
        return smiles_from_json
    else:
        raise ValueError("Invalid input")
    


# In[ ]:





# In[11]:


def cid_to_SMILES(input_id: int)-> str:
    """
    Convert PubChem id into SMILES notation.

    Parameters:
    input_id (str): The PubChem id of the chemical compound.

    Returns:
    str: The SMILES notation in the output format.
    """
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
    """
    Convert molecular formula into SMILES notation.

    Parameters:
    input_formula (str): The molecular formula of the chemical compound.

    Returns:
    str: The SMILES notation in the output format.
    """
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
    """
    Convert zinc id into SMILES notation.

    Parameters:
    input_id (str): The ZINC15 id of the chemical compound.

    Returns:
    str: The SMILES notation in the output format.
    """   
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



