from datetime import timedelta

from fastapi import APIRouter, Depends, HTTPException, status
from fastapi.security import OAuth2PasswordRequestForm
from sqlalchemy.ext.asyncio import AsyncSession

from app.api.dependencies.database import get_async_db
from app.core.security import create_access_token, verify_password
from app.schemas import token as token_schema
from app.schemas import users as user_schema
from app.services import users as user_service
from app.settings import settings

router = APIRouter()


@router.post("/register", response_model=user_schema.User)
async def register_user(
    user_in: user_schema.UserCreate, db: AsyncSession = Depends(get_async_db)
):
    """
    Create a new user.
    """
    user = await user_service.get_user_by_email(db, email=user_in.email)
    if user:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="The user with this email is already registered.",
        )
    user = await user_service.create_user(db, user_in=user_in)
    return user


@router.post("/login/access-token", response_model=token_schema.Token)
async def login_for_access_token(
    db: AsyncSession = Depends(get_async_db),
    form_data: OAuth2PasswordRequestForm = Depends(),
):
    """
    OAuth2 compatible token login, get an access token for future requests.
    """
    user = await user_service.get_user_by_email(db, email=form_data.username)
    if not user or not verify_password(form_data.password, user.hashed_password):
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Incorrect email or password",
            headers={"WWW-Authenticate": "Bearer"},
        )
    access_token_expires = timedelta(minutes=settings.ACCESS_TOKEN_EXPIRE_MINUTES)
    access_token = create_access_token(
        user.id, expires_delta=access_token_expires
    )
    return {
        "access_token": access_token,
        "token_type": "bearer",
    }
