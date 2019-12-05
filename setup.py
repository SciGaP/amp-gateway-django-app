import setuptools

setuptools.setup(
    # Technically, not a django app, yet. Only has output view providers.
    name="amp-gateway-django-app",
    version="0.0.1",
    description="Custom output viewer for AMP Gateway",
    packages=setuptools.find_packages(),
    install_requires=[
        'numpy',
        'matplotlib',
        'scipy'
    ],
    entry_points="""
[airavata.output_view_providers]
trecx-plot = amp_gateway.plot:TRecXPlotViewProvider
""",
)
