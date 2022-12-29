#!python
import os

def run(port=10077):
    print('Starting PeptDeep Web Server ...')

    _this_file = __file__
    _this_directory = os.path.dirname(_this_file)

    file_path = os.path.join(_this_directory, 'webui', 'main_ui.py')

    HOME = os.path.expanduser("~")

    ST_PATH = os.path.join(HOME, ".streamlit")

    if not os.path.isdir(ST_PATH):
        os.mkdir(ST_PATH)

    #Check if streamlit credentials exists
    ST_CREDENTIALS = os.path.join(ST_PATH, 'credentials.toml')
    if not os.path.isfile(ST_CREDENTIALS):
        with open(ST_CREDENTIALS, 'w') as file:
            file.write("[general]\n")
            file.write('\nemail = ""')


    import sys
    from streamlit.web import cli as stcli

    theme = []

    theme.append("--theme.backgroundColor=#FFFFFF")
    theme.append("--theme.secondaryBackgroundColor=#f0f2f6")
    theme.append("--theme.textColor=#262730")
    theme.append("--theme.font=sans serif")
    theme.append("--theme.primaryColor=#18212b")

    args = [
        "streamlit", "run", 
        file_path, "--global.developmentMode=false", 
        f"--server.port={port}", 
        "--browser.gatherUsageStats=False",
        "--logger.level=error"
    ]

    args.extend(theme)

    sys.argv = args

    sys.exit(stcli.main())