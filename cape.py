import streamlit as st
from io import StringIO
from Bio.PDB import PDBParser
import pandas as pd
import requests
import py3Dmol
from chembl_webresource_client.new_client import new_client

st.title("Protein-Ligand Insight")
st.write("A tool for exploring protein-ligand interactions, PubChem compounds, and FDA-approved drugs in ChEMBL.")

# =============================================================================
# Protein Neighbors Functionality
# =============================================================================
st.header("Protein Neighbors")
st.write("Find coordinates of residues around ligands/ions in a protein PDB file.")

protein_pdb_file = st.file_uploader("Upload protein PDB file", type="pdb", key="protein_file_uploader")
distance = st.number_input("Distance cutoff (Å)", min_value=0.1, value=5.0, key="protein_distance_input")

def get_coordinates_around(structure, distance=5.0):
    close_residues = []
    ligand_or_ion_residues = []
    ligand_names = {}

    for model in structure:
        for chain in model.get_chains():
            for residue in chain.get_residues():
                if residue.id[0] != ' ':
                    ligand_or_ion_residues.append(residue)
                    ligand_names[residue.id[1]] = residue.resname

    if not ligand_or_ion_residues:
        st.warning("No ligands or ions found in the structure.")
        return pd.DataFrame()

    for model in structure:
        for chain in model.get_chains():
            for residue in chain.get_residues():
                if residue.id[0] == ' ':
                    for ligand_residue in ligand_or_ion_residues:
                        for atom in residue.get_atoms():
                            for ligand_atom in ligand_residue.get_atoms():
                                if atom - ligand_atom < distance:
                                    close_residues.append((residue, ligand_residue.id[1]))
                                    break
                            else:
                                continue
                            break
                    else:
                        continue
                    break

    close_residues = list(set(close_residues))

    data = []
    for residue, ligand_id in close_residues:
        for atom in residue.get_atoms():
            data.append({
                'chain_id': residue.parent.id,
                'residue_number': residue.id[1],
                'residue_name': residue.resname,
                'atom_name': atom.name,
                'x': atom.coord[0],
                'y': atom.coord[1],
                'z': atom.coord[2],
                'ligand_ion_name': ligand_names.get(ligand_id, "Unknown"),
            })
    return pd.DataFrame(data)

if protein_pdb_file:
    protein_pdb_file_content = protein_pdb_file.read().decode("utf-8")
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("protein", StringIO(protein_pdb_file_content))
        df = get_coordinates_around(structure, distance)
        if not df.empty:
            st.subheader(f"Residue coordinates within {distance} Å of ligand/ions:")
            st.dataframe(df)
        else:
            st.warning(f"No ligands or ions found, or no residues within {distance} Å.")
    except Exception as e:
        st.error(f"Error parsing PDB file: {e}")

# =============================================================================
# PubChem Functionality
# =============================================================================
st.header("PubChem Compound Explorer")
st.write("Search PubChem for compounds, view their 3D structures, and download SDF files.")

search_query = st.text_input("Enter compound name or PubChem CID:", key="pubchem_search_input")

def search_pubchem(query):
    base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
    output = "cids/JSON"
    url = f"{base_url}{query}/{output}"
    try:
        response = requests.get(url)
        response.raise_for_status()
        data = response.json()
        if 'IdentifierList' in data and 'CID' in data['IdentifierList']:
            return data['IdentifierList']['CID']
        else:
            return []
    except requests.exceptions.RequestException as e:
        st.error(f"Error searching PubChem: {e}")
        return []

def get_sdf_3d(cid):
    base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
    output = "SDF"
    record_type = "record_type=3d"
    url = f"{base_url}{cid}/{output}?{record_type}"
    try:
        response = requests.get(url)
        response.raise_for_status()
        return response.content
    except requests.exceptions.RequestException as e:
        st.error(f"Error retrieving 3D SDF for CID {cid}: {e}")
        return None

def display_molecule_3d(sdf_content):
    sdf_str = sdf_content.decode('utf-8')
    viewer = py3Dmol.view(width=400, height=400)
    viewer.addModel(sdf_str, 'sdf')
    viewer.setStyle({'stick': {}})
    viewer.zoomTo()
    html = viewer._repr_html_()
    st.components.v1.html(html, height=400)

if st.button("Search", key="pubchem_search_button"):
    if search_query:
        cids = search_pubchem(search_query)
        if cids:
            st.subheader(f"Found {len(cids)} compounds:")
            selected_cid = st.selectbox("Select a CID to view 3D structure and download SDF:", cids, key="pubchem_cid_select")
            if selected_cid:
                sdf_content = get_sdf_3d(selected_cid)
                if sdf_content:
                    st.subheader(f"3D Structure of CID: {selected_cid}")
                    display_molecule_3d(sdf_content)
                    st.download_button(
                        label="Download 3D SDF File",
                        data=sdf_content,
                        file_name=f"{selected_cid}.sdf",
                        mime="chemical/x-mdl-sdfile"
                    )
                else:
                    st.warning(f"No 3D structure available for CID: {selected_cid}")
        else:
            st.warning(f"No compounds found matching '{search_query}'")
    else:
        st.warning("Please enter a compound name or CID to search.")

# =============================================================================
# ChEMBL Drug Search Section
# =============================================================================
st.header("FDA Approved Drug Search (ChEMBL)")
drug_name = st.text_input("Enter drug name:", key="chembl_drug_name_input")

def search_chembl(drug_name):
    molecule = new_client.molecule
    results = molecule.filter(max_phase=4, pref_name__icontains=drug_name)
    return results

if st.button("Search", key="chembl_search_button"):
    if drug_name:
        st.info(f"Searching ChEMBL for FDA approved drugs containing '{drug_name}'...")
        results = list(search_chembl(drug_name))
        if results:
            st.success(f"Found {len(results)} FDA approved drug(s) containing '{drug_name}':")
            for r in results:
                st.subheader(f"ChEMBL ID: {r['molecule_chembl_id']}")
                st.write(f"**Preferred Name:** {r['pref_name']}")
                synonyms = r.get('molecule_synonyms')
                if synonyms:
                    string_synonyms = [syn['synonym'] for syn in synonyms if isinstance(syn, dict) and 'synonym' in syn and isinstance(syn['synonym'], str)]
                    st.write(f"**Synonyms:** {', '.join(string_synonyms) if string_synonyms else 'N/A'}")
                else:
                    st.write("**Synonyms:** N/A")
                st.write(f"**Highest Phase:** {r['max_phase']}")
                st.write("-" * 20)
        else:
            st.warning(f"No FDA approved drugs found containing '{drug_name}'.")
    else:
        st.warning("Please enter a drug name to search.")