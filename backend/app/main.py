import os

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from app.db import init_db
from app.routers import compound, health, predict, training


app = FastAPI(
    title="MRGPRX2 Ligand Predictor API",
    version="0.1.0",
    description="Backend API for compound lookup and MRGPRX2 prediction workflows.",
)


@app.on_event("startup")
def _create_tables() -> None:
    init_db()

def _parse_origins(raw: str) -> list[str]:
    """Splits a comma-separated FRONTEND_ORIGIN value into trimmed,
    slash-stripped origins - so pasted values with stray whitespace or a
    trailing "/" (a common copy-paste mistake with deployment URLs) still
    match, and multiple origins (e.g. a Vercel production domain plus a
    preview deployment URL) can be allow-listed at once.
    """
    return [origin.strip().rstrip("/") for origin in raw.split(",") if origin.strip()]


frontend_origins = _parse_origins(os.getenv("FRONTEND_ORIGIN", "http://localhost:3000"))

app.add_middleware(
    CORSMiddleware,
    allow_origins=frontend_origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

app.include_router(health.router)
app.include_router(compound.router)
app.include_router(predict.router)
app.include_router(training.router)
