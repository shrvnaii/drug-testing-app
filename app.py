import streamlit as st
import urllib.request
import ssl
import py3Dmol
from stmol import showmol

st.title("In Silico Drug Tester")
st.write("Analyze small molecules for bioavailability and visualize 3D structures.")

smiles_input = st.text_input("Enter Molecule SMILES string:", "CC(=O)OC1=CC=CC=C1C(=O)O")

if smiles_input:
    mol = Chem.MolFromSmiles(smiles_input)
    
    if mol:
        st.subheader("Lipinski's Rule of 5 (Bioavailability)")
        
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        
        col1, col2, col3, col4 = st.columns(4)
        col1.metric("Weight (< 500 Da)", f"{mw:.2f}")
        col2.metric("LogP (< 5)", f"{logp:.2f}")
        col3.metric("H-Donors (< 5)", f"{hbd}")
        col4.metric("H-Acceptors (< 10)", f"{hba}")
        
        st.subheader("3D Structure Viewer")
        # --- NEW SCREENING MODULE START ---
        st.subheader("Safety & Environmental Screening")
        
        # Define SMARTS patterns for structural alerts
        alerts = {
            "Halogenated Groups (Low Biodegradability)": "[F,Cl,Br,I]",
            "Epoxides (Skin Sensitization/Highly Reactive)": "C1OC1",
            "Heavy Metals (High Toxicity)": "[As,Sb,Hg,Pb,Cd]",
            "Phenols (Potential Irritant)": "c1ccccc1O",
            "Aldehydes (Skin Sensitization)": "[CX3H1](=O)[#6]"
        }
        
        flags = []
        for alert_name, smarts in alerts.items():
            pattern = Chem.MolFromSmarts(smarts)
            if pattern and mol.HasSubstructMatch(pattern):
                flags.append(alert_name)
                
        if flags:
            st.warning("⚠️ Structural Alerts Detected:")
            for flag in flags:
                st.write(f"- {flag}")
            st.info("Note: Presence of these groups requires further literature review for safe usage.")
        else:
            st.success("✅ No common toxic or non-biodegradable structural alerts detected.")
        st.divider()
        # --- NEW SCREENING MODULE END ---
        
        mol_3d = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol_3d)
        AllChem.MMFFOptimizeMolecule(mol_3d)
        
        mol_block = Chem.MolToMolBlock(mol_3d)
        
        view = py3Dmol.view(width=700, height=500)
        view.addModel(mol_block, "mol")
        view.setStyle({"stick": {}, "sphere": {"scale": 0.3}})
        view.zoomTo()
        
        showmol(view, height=500, width=700)
    else:
        st.error("Invalid SMILES string. Please check your input.")
# --- NEW ENZYME VIEWER MODULE START ---
st.divider()
st.subheader("Industrial Enzyme & Protein Viewer")
st.write("Enter a 4-letter Protein Data Bank (PDB) ID to visualize a macromolecule.")

# Default is 1BAG (an Alpha-Amylase)
pdb_id = st.text_input("Enter PDB ID:", "1BAG")

# NEW: Toggle switch for the surface map
show_surface = st.checkbox("Show Molecular Surface Map (Identifies Binding Pockets)")

if pdb_id:
    try:
        # Fetch the protein structure directly from the RCSB PDB database
        url = f"https://files.rcsb.org/view/{pdb_id.upper()}.pdb"
        ssl_context = ssl._create_unverified_context()
        req = urllib.request.Request(url)
        response = urllib.request.urlopen(req, context=ssl_context)
        pdb_data = response.read().decode('utf-8')
        
        # Set up a new 3D viewer for the large protein
        view2 = py3Dmol.view(width=700, height=500)
        view2.addModel(pdb_data, "pdb")
        
        # Style the protein as a cartoon ribbon, colored by structure
        view2.setStyle({'protein': {}}, {'cartoon': {'color': 'spectrum'}})
        
        # If there is a small substrate/ligand bound to the enzyme, highlight it in sticks
        view2.setStyle({'ligand': {}}, {'stick': {'colorscheme': 'greenCarbon', 'radius': 0.2}})
        
        # NEW: Calculate and render the van der Waals surface if checkbox is ticked
        if show_surface:
            view2.addSurface(py3Dmol.VDW, {'opacity': 0.6, 'color': 'lightblue'}, {'protein': {}})
        
        view2.zoomTo()
        showmol(view2, height=500, width=700)
        
    except Exception as e:
        st.error("Could not fetch that PDB ID. Make sure it is a valid 4-letter code.")
# --- NEW ENZYME VIEWER MODULE END ---
