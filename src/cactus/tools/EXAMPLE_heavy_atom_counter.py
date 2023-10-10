import pubchempy
from langchain.agents import Tool

def heavy_atom_counter(compound:str) -> int:
    '''
    mockup function: given a compound name,
    return heavy atom count as an integer.

    parameters:
        compound: a chemical compound name (str)

    returns:
        heavy_atom_count (int)

    '''

    Query_Compound = pubchempy.get_compounds(str(compound),'name')[0]
    Heavy_Atom_Count = Query_Compound.heavy_atom_count

    return Heavy_Atom_Count

Heavy_Atom = Tool(
	name="heavy_atom_counter",
	func=heavy_atom_counter,
	description="Helps to retrieve the number of heavy atoms present in a compound. Input should be the name of a chemical ('compound') as a string.
)
