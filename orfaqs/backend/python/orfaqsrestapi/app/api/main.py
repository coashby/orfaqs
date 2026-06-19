"""
Main application module for the ORFaqs REST API.

This module initializes and configures the FastAPI application instance,
including route registration and API prefix setup for the ORFaqs REST API.
"""

from fastapi import FastAPI

from orfaqs.backend.python.orfaqsrestapi.app.api.configuration import (
    orfaqs_rest_api_settings,
)
from orfaqs.backend.python.orfaqsrestapi.app.api.routes import discoverproteins

app = FastAPI()
app.include_router(
    discoverproteins.router,
    prefix=orfaqs_rest_api_settings.ROUTE_PREFIX,
)
