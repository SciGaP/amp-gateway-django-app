# AMP Gateway Django Portal Integrations

## Getting Started

1. Install this Python package into the Airavata Django Portal virtual
   environment.

   ```
   pip install git+https://github.com/SciGaP/amp-gateway-django-app.git#egg=amp_gateway_django_app
   ```

2. Configure application output file's metadata to have the `trecx-plot` output
   view provider. Under the **Application Interface**, configure the
   **Metadata** of the output field to have the following:
   ```json
   {
     "output-view-providers": ["trecx-plot"]
   }
   ```

## Development

See https://github.com/apache/airavata-django-portal for instructions on setting
up a local Airavata Django Portal.

For local development, activate the Airavata Django Portal virtual environment,
then in this project's root directory, run the following:

```
pip install -r requirements.txt
python setup.py develop
```
