#copy from alphapept
import streamlit as st
from peptdeep.webui.ui_utils import markdown_link
import peptdeep


def show():
    """Streamlit page that displays information on how to get started."""
    st.write("# Start page")
    st.write("Welcome to AlphaPeptDeep (PeptDeep for short).")
    st.write("#### Cite: ")
    st.write("Wen-Feng Zeng, Xie-Xuan Zhou, Sander Willems, Constantin Ammar, Maria Wahle, Isabell Bludau, Eugenia Voytik, Maximillian T. Strauss & Matthias Mann. AlphaPeptDeep: a modular deep learning framework to predict peptide properties for proteomics. Nat Commun 13, 7238 (2022). https://doi.org/10.1038/s41467-022-34904-3")

    with st.expander("Navigation"):
        st.write("Use the sidebar on the left to navigate through the different tasks.")
        st.write(
            "- Model: Set the settings of peptdeep models."
            " \n- Transfer: Transfer learning and save the refined models."
            " \n- Library: Predict spectral libraries."
            " \n- Rescore: Rescore DDA search results (the GUI is not finished yet)."
            " \n- Server: Check the job queue for `Transfer`, `Library`, ... tasks."
            " \n- Settings: Save settings or upload previous settings in YAML files."
        )

    with st.expander("Resources"):
        st.write(
            "On the following pages you can find additional information about PeptDeep:"
        )
        markdown_link("PeptDeep on GitHub", peptdeep.__github__)
        markdown_link("Code documentation", peptdeep.__doc__)
        markdown_link(
            "Report an issue or requirement a feature (GitHub)",
            f"{peptdeep.__github__}/issues/new/choose",
        )
        markdown_link(
            "Contact (e-mail)",
            f"mailto:{peptdeep.__author_email__}?subject=peptdeep (v{peptdeep.__version__})",
        )

    with st.expander("Server"):

        st.write(
            "When starting PeptDeep you launch a server that can be closed when closing the terminal window."
        )
        st.write(
            "If your firewall policy allows external access, this page can also be accessed from other computers in the network."
        )
        st.write(
            "The server starts an PeptDeep process in the background that will process new experiments once they are submitted."
        )