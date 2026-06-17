"""Database engine/session setup for persisted application data.

Reads `DATABASE_URL` so the same code runs against a local SQLite file in
development and a hosted Postgres instance (e.g. Supabase/Neon) in
production - hosted Postgres survives redeploys on PaaS platforms whose
filesystems are wiped on each deploy, whereas a SQLite file on disk would
not. Defaults to a SQLite file under `backend/data/` when unset.
"""

from __future__ import annotations

import os
from pathlib import Path

from dotenv import load_dotenv
from sqlalchemy import create_engine
from sqlalchemy.orm import DeclarativeBase, Session, sessionmaker

load_dotenv(Path(__file__).resolve().parent.parent / ".env")

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
DEFAULT_DATABASE_URL = f"sqlite:///{DATA_DIR / 'app.db'}"

DATABASE_URL = os.getenv("DATABASE_URL", DEFAULT_DATABASE_URL)

_connect_args = {"check_same_thread": False} if DATABASE_URL.startswith("sqlite") else {}

if DATABASE_URL.startswith("sqlite"):
    DATA_DIR.mkdir(parents=True, exist_ok=True)

try:
    engine = create_engine(DATABASE_URL, connect_args=_connect_args)
except Exception as _e:  # pragma: no cover
    import logging
    logging.getLogger(__name__).warning(
        "DATABASE_URL invalid (%s); falling back to SQLite at %s",
        _e, DEFAULT_DATABASE_URL,
    )
    engine = create_engine(DEFAULT_DATABASE_URL, connect_args={"check_same_thread": False})

SessionLocal = sessionmaker(bind=engine, autoflush=False, expire_on_commit=False)


class Base(DeclarativeBase):
    pass


def init_db() -> None:
    from app import models  # noqa: F401  (ensure models are registered before create_all)

    Base.metadata.create_all(bind=engine)


def get_session() -> Session:
    return SessionLocal()
