from sqlalchemy import Column, Integer, String
from .database import Base


class Molecule(Base):
    __tablename__ = "molecules"

    id = Column(Integer, primary_key=True, index=True)
    structure = Column(String, nullable=False)
