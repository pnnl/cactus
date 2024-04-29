import pubchempy as pcp


def convert_identifiers(input_name: str, input_type: str, output_type: str) -> str:
    """
    Convert chemical identifiers between different types.

    Parameters:
    input_name (str): The identifier of the chemical.
    input_type (str): The type of the identifier being input.
    output_type (str): The type of the identifier to output.

    Returns:
    str: The identifier in the output format.
    """

    # Instantiate compound object directly if input_type is 'CID'
    if input_type.lower() == "cid":
        compound = pcp.Compound.from_cid(int(input_name))
    else:
        # Get the first compound that matches the input_name and input_type
        compound = pcp.get_compounds(input_name, input_type)[0]

    # Get the first common synonym for the name if output_type is 'name'
    if output_type.lower() == "name":
        return compound.synonyms[0] if compound.synonyms else "No common name found"

    # Otherwise, get the attribute specified by output_type
    return getattr(compound, output_type, "Invalid output type")


# Example usage
smiles = convert_identifiers("Aspirin", "name", "isomeric_smiles")
print(smiles)
