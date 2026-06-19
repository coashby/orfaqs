"""
Configuration module for the ORFaqs REST API.

This module defines settings and configuration for the ORFaqs REST API
application, including API versioning, route naming, and environment variable
management using Pydantic settings.
"""

from pydantic_settings import BaseSettings, SettingsConfigDict


class ORFaqsRestApiSettings(BaseSettings):
    model_config = SettingsConfigDict(
        env_file='../.env',
        env_ignore_empty=True,
        extra='ignore',
    )
    API_V1_STR: str = '/api/v1'
    ROUTE_NAME: str = 'orfaqs'
    ROUTE_PREFIX: str = f'/{ROUTE_NAME}'


orfaqs_rest_api_settings = ORFaqsRestApiSettings()
