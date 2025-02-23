from sqlalchemy import Column, Integer, String

try:
    from src.database import Base
except ModuleNotFoundError:
    from database import Base


class Molecule(Base):
    __tablename__ = "molecules"

    id = Column(Integer, primary_key=True, index=True)
    structure = Column(String, nullable=False)
