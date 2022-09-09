#copy from alphapept
import streamlit as st
from peptdeep.webui.ui_utils import markdown_link
import peptdeep


def show():
    """Streamlit page that displays information on how to get started."""
    st.write("# Getting started")
    st.text("Welcome to AlphaPeptDeep (PeptDeep for short).")

    with st.expander("Navigation"):
        st.write("Use the sidebar on the left to navigate through the different menus.")
        st.write(
            "- Server: Start/stop the server process to run `transfer learning`, `library prediction`, ... tasks."
            " \n- Settings: Save settings or upload previous settings."
            " \n- Model: Set the settings of peptdeep models."
            " \n- Transfer: Transfer learning and save the refined models."
            " \n- Library: Predict spectral libraries."
            " \n- Rescore: Rescore DDA search results (the GUI is not finished yet)."
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

    with st.expander("Sample Run"):
        st.write("To Be Finished")