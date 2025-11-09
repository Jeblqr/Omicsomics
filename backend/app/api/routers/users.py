"""User-related API endpoints."""

from fastapi import APIRouter, Depends

from app.api.dependencies.auth import get_current_user
from app.models.user import User
from app.schemas import users as user_schema

router = APIRouter(tags=["users"])


@router.get("/me", response_model=user_schema.User)
async def read_current_user(current_user: User = Depends(get_current_user)) -> user_schema.User:
    """Return the currently authenticated user."""
    return current_user
