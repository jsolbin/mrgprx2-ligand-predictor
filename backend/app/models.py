from __future__ import annotations

from datetime import datetime, UTC

from sqlalchemy import DateTime, String, Text
from sqlalchemy.orm import Mapped, mapped_column

from app.db import Base


class LabeledCompound(Base):
    __tablename__ = "labeled_compounds"

    id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    name: Mapped[str] = mapped_column(String(255), nullable=False)
    smiles: Mapped[str] = mapped_column(String(1024), nullable=False, unique=True, index=True)
    label: Mapped[str] = mapped_column(String(32), nullable=False)
    note: Mapped[str | None] = mapped_column(Text, nullable=True)
    source: Mapped[str] = mapped_column(String(32), nullable=False, default="user")
    submitted_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), nullable=False, default=lambda: datetime.now(UTC)
    )

    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "smiles": self.smiles,
            "label": self.label,
            "note": self.note,
            "source": self.source,
            "submitted_at": self.submitted_at.isoformat(),
        }
