import os

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from app.routers import compound, health, predict, training


app = FastAPI(
    title="MRGPRX2 Ligand Predictor API",
    version="0.1.0",
    description="Backend API for compound lookup and MRGPRX2 prediction workflows.",
)

frontend_origin = os.getenv("FRONTEND_ORIGIN", "http://localhost:3000")

app.add_middleware(
    CORSMiddleware,
    allow_origins=[frontend_origin],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

app.include_router(health.router)
app.include_router(compound.router)
app.include_router(predict.router)
app.include_router(training.router)
