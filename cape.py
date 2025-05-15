import streamlit as st
from io import StringIO
from Bio.PDB import PDBParser
import pandas as pd
import requests
import py3Dmol
from chembl_webresource_client.new_client import new_client
import plotly.express as px

st.title("Protein-Ligand Insight")
st.write(
    "A tool for exploring protein-ligand interactions, PubChem compounds, and FDA-approved drugs in ChEMBL."
)

# =============================================================================
# Protein Neighbors Functionality
# =============================================================================
st.header("Protein Neighbors")
st.write("Find coordinates of residues around ligands/ions in a protein PDB file.")

protein_pdb_file = st.file_uploader(
    "Upload protein PDB file", type="pdb", key="protein_file_uploader"
)
distance = st.number_input(
    "Distance cutoff (Å)", min_value=0.1, value=5.0, key="protein_distance_input"
)


def get_coordinates_around(structure, distance=5.0):
    close_residues = []
    ligand_or_ion_residues = []
    ligand_names = {}

    for model in structure:
        for chain in model.get_chains():
            for residue in chain.get_residues():
                if residue.id[0] != " ":
                    ligand_or_ion_residues.append(residue)
                    ligand_names[residue.id[1]] = residue.resname

    if not ligand_or_ion_residues:
        st.warning("No ligands or ions found in the structure.")
        return pd.DataFrame()

    for model in structure:
        for chain in model.get_chains():
            for residue in chain.get_residues():
                if residue.id[0] == " ":
                    for ligand_residue in ligand_or_ion_residues:
                        for atom in residue.get_atoms():
                            for ligand_atom in ligand_residue.get_atoms():
                                if atom - ligand_atom < distance:
                                    close_residues.append(
                                        (residue, ligand_residue.id[1])
                                    )
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
            data.append(
                {
                    "chain_id": residue.parent.id,
                    "residue_number": residue.id[1],
                    "residue_name": residue.resname,
                    "atom_name": atom.name,
                    "x": atom.coord[0],
                    "y": atom.coord[1],
                    "z": atom.coord[2],
                    "ligand_ion_name": ligand_names.get(ligand_id, "Unknown"),
                }
            )
    return pd.DataFrame(data)


if protein_pdb_file:
    protein_pdb_file_content = protein_pdb_file.read().decode("utf-8")
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("protein", StringIO(protein_pdb_file_content))
        df = get_coordinates_around(structure, distance)
        if not df.empty:
            st.subheader(
                f"Residue coordinates within {distance} Å of ligand/ions:"
            )

            # Visualization (Example - 3D Scatter Plot)
            fig = px.scatter_3d(
                df,
                x="x",
                y="y",
                z="z",
                color="residue_name",
                hover_data=["chain_id", "residue_number", "atom_name", "ligand_ion_name"],
                title="3D View of Residues Around Ligand/Ion",
            )
            st.plotly_chart(fig)

            # Display the table
            st.dataframe(df)
        else:
            st.warning(
                f"No ligands or ions found, or no residues within {distance} Å."
            )
    except Exception as e:
        st.error(f"Error parsing PDB file: {e}")

# =============================================================================
# PubChem Functionality
# =============================================================================
st.header("PubChem Compound Explorer")
st.write(
    "Search PubChem for compounds, view their 3D structures, and download SDF files."
)

search_query = st.text_input(
    "Enter compound name or PubChem CID:", key="pubchem_search_input"
)


@st.cache_data
def search_pubchem(query):
    base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
    output = "cids/JSON"
    url = f"{base_url}{query}/{output}"
    try:
        response = requests.get(url)
        response.raise_for_status()
        data = response.json()
        if "IdentifierList" in data and "CID" in data["IdentifierList"]:
            return data["IdentifierList"]["CID"]
        else:
            return []
    except requests.exceptions.RequestException as e:
        st.error(f"Error searching PubChem: {e}")
        return []


@st.cache_resource
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
    sdf_str = sdf_content.decode("utf-8")
    viewer = py3Dmol.view(width=400, height=400)
    viewer.addModel(sdf_str, "sdf")
    viewer.setStyle({"stick": {}})
    viewer.zoomTo()
    html = viewer.render_html().strip()
    st.components.v1.html(f"<div>{html}</div>", height=400)


if st.button("Search", key="pubchem_search_button"):
    if search_query:
        cids = search_pubchem(search_query)
        if cids:
            st.subheader(f"Found {len(cids)} compounds:")
            selected_cid = st.selectbox(
                "Select a CID to view 3D structure and download SDF:",
                cids,
                key="pubchem_cid_select",
            )
            if selected_cid:
                sdf_content = get_sdf_3d(selected_cid)
                if sdf_content:
                    st.subheader(f"3D Structure of CID: {selected_cid}")
                    display_molecule_3d(sdf_content)
                    st.download_button(
                        label="Download 3D SDF File",
                        data=sdf_content,
                        file_name=f"{selected_cid}.sdf",
                        mime="chemical/x-mdl-sdfile",
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


@st.cache_data
def search_chembl(drug_name):
    molecule = new_client.molecule
    results = molecule.filter(max_phase=4, pref_name__icontains=drug_name).only(
        ["molecule_chembl_id", "pref_name", "max_phase", "therapeutic_flag", "molecule_type"]
    )
    return pd.DataFrame(results)


if st.button("Search", key="chembl_search_button"):
    if drug_name:
        st.info(f"Searching ChEMBL for FDA approved drugs containing '{drug_name}'...")
        results_df = search_chembl(drug_name)
        if not results_df.empty:
            st.success(
                f"Found {len(results_df)} FDA approved drug(s) containing '{drug_name}':"
            )

            st.subheader("Distribution of Molecule Types")
            fig_type = px.bar(
                results_df,
                x="molecule_type",
                title="Molecule Types of Found Drugs",
            )
            st.plotly_chart(fig_type)

            st.subheader("Distribution of Drug Development Phases")
            fig_phase = px.histogram(
                results_df,
                x="max_phase",
                title="Drug Development Phases (All are Phase 4)",
            )
            st.plotly_chart(fig_phase)

            st.subheader("Drug Information")
            st.dataframe(results_df)

            for index, row in results_df.iterrows():
                st.markdown(f"**ChEMBL ID:** {row['molecule_chembl_id']}")
                st.markdown(f"**Preferred Name:** {row['pref_name']}")
                st.markdown(f"**Highest Phase:** {row['max_phase']}")
                st.markdown(f"**Therapeutic Flag:** {row['therapeutic_flag']}")
                st.markdown(f"**Molecule Type:** {row['molecule_type']}")
                st.markdown("---")
        else:
            st.warning(f"No FDA approved drugs found containing '{drug_name}'.")
    else:
        st.warning("Please enter a drug name to search.")

elif section == "About":
    st.header("About Protein-Ligand Insight")
    st.markdown(
        """
        Welcome to Protein-Ligand Insight! This application provides a suite of tools
        for researchers and enthusiasts interested in the world of biomolecules and drug discovery.

        **Our Mission:**
        To provide an intuitive and informative platform for exploring the interactions
        between proteins and small molecules, as well as accessing valuable information
        about chemical compounds and approved drugs.

        **Key Features:**
        """
    )
    
    features = [
        "Protein Neighbors: Upload your own Protein Data Bank (PDB) files. Specify a distance cutoff to identify amino acid residues in the vicinity of ligands or ions within the protein structure. Visualize the spatial arrangement of these neighboring residues in 3D. Download the list of neighboring residues and their coordinates.",
        "PubChem Explorer: Search the vast PubChem database by compound name or Chemical Identifier (CID). View interactive 3D structures of retrieved compounds. Download the 3D structural data in SDF format.",
        "ChEMBL Drug Search: Explore FDA-approved drugs using the ChEMBL database. Search for drugs by name and view their development phase and other relevant information. Visualize the types of molecules found in your search.",
    ]
    
    for feature in features:
        st.markdown(f"- {feature}")

    st.markdown(
        """

        **Technology Stack:**
        - Streamlit: For creating the interactive web application.
        - Biopython: For parsing and analyzing PDB files.
        - Pandas: For data manipulation and analysis.
        - Requests: For making HTTP requests to external APIs (PubChem).
        - py3Dmol: For rendering 3D molecular structures.
        - ChEMBL Web Resource Client: For interacting with the ChEMBL database.
        - Plotly Express: For creating interactive visualizations.

        **Acknowledgements:**

        This project's successful completion owes much to the collaborative spirit and support I received.
        I am deeply indebted to Dr. Kushagra Kashyap, for his invaluable mentorship, constructive feedback, and ongoing motivation.
        His expertise in bioinformatics was essential to the development of this study.
        I also acknowledge and appreciate the stimulating academic environment and helpful discussions
        provided by the faculty and my fellow students in the Department of Bioinformatics.

        **Developed By:** [Yash Birjajdar/DESPU]
        **Contact:** [3522411027@despu.edu.in/Link]
        **Version:** [Streamlit v1.45.1]
        """
    )
    st.markdown("---")
    st.markdown(
        """
        **GitHub Profile:** [https://github.com/YASHBIOINFOWORK/bio_info.git]
        """
    )
