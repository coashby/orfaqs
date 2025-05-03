# Setting Up Your Dev Environment
## Recommended Components
- Python 3.12
- Multi-core CPU

## Basic Setup
1. Clone the repo and navigate to the `/python` directory.
   ```shell
   git clone https://github.com/coashby/orfaqs.git && \
   cd orfaqs/python/
   ```
1. Open your preferred code editor in this directory.
1. Create a virtual environment.
    ```shell
    python3 -m venv .venv && source .venv/bin/activate
    ```
1. Add the current directory to your python path environment variable.
    ```
    export PYTHONPATH=$(pwd)
    ```
1. Make sure `pip` is up-to-date and install the necessary library dependencies
using the requirements file provided.
    ```shell
     python -m pip install --upgrade pip && pip install -r ./sop/requirements.txt
    ````