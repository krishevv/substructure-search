from os import getenv
from fastapi import FastAPI, Depends, HTTPException, UploadFile, File
from sqlalchemy.orm import Session
from rdkit import Chem
import json
import os

from .database import get_db, Base, engine
from .models import Molecule
from pydantic import BaseModel


# Инициализация приложения FastAPI
app = FastAPI()

# Создание таблиц, если они не существуют

if os.getenv("RUN_MIGRATIONS", "false").lower() == "true":
    Base.metadata.create_all(bind=engine)


class MoleculeCreate(BaseModel):
    structure: str


class MoleculeUpdate(BaseModel):
    structure: str


class SubstructureQuery(BaseModel):
    substructure: str

# Эндпоинты


@app.get("/")
def get_server():
    return {"server_id": getenv("SERVER_ID", "unknown")}


@app.post("/add")
def add_molecule(molecule: MoleculeCreate, db: Session = Depends(get_db)):
    db_molecule = Molecule(structure=molecule.structure)
    db.add(db_molecule)
    db.commit()
    db.refresh(db_molecule)
    return {"message": "Молекула добавлена", "id": db_molecule.id}


@app.get("/molecule/{id}")
def get_molecule(id: int, db: Session = Depends(get_db)):
    molecule = db.query(Molecule).filter(Molecule.id == id).first()
    if not molecule:
        raise HTTPException(status_code=404, detail="Молекула не найдена")
    return molecule


@app.put("/molecule/{id}")
def update_molecule(id: int, molecule: MoleculeUpdate, db: Session = Depends(get_db)):
    db_molecule = db.query(Molecule).filter(Molecule.id == id).first()
    if not db_molecule:
        raise HTTPException(status_code=404, detail="Молекула не найдена.")
    db_molecule.structure = molecule.structure
    db.commit()
    db.refresh(db_molecule)

    return {"message": "Молекула обновлена успешно.", "id": db_molecule.id}


@app.delete("/molecule/{id}")
def delete_molecule(id: int, db: Session = Depends(get_db)):
    molecule = db.query(Molecule).filter(Molecule.id == id).first()
    if not molecule:
        raise HTTPException(status_code=404, detail="Молекула не найдена")
    db.delete(molecule)
    db.commit()
    return {"message": "Молекула удалена", "id": id}


@app.get("/list")
def list_molecules(db: Session = Depends(get_db)):
    return db.query(Molecule).all()


@app.post("/search")
def search_substructure(query: SubstructureQuery, db: Session = Depends(get_db)):
    substructure = Chem.MolFromSmiles(query.substructure)
    if not substructure:
        raise HTTPException(status_code=400, detail="Некорректный SMILES")

    matching_molecules = []
    for molecule in db.query(Molecule).all():
        mol = Chem.MolFromSmiles(molecule.structure)
        if mol and mol.HasSubstructMatch(substructure):
            matching_molecules.append({"id": molecule.id, "smiles": molecule.structure})

    return {"matches": matching_molecules}


@app.post("/upload")
async def upload_file(file: UploadFile = File(...), db: Session = Depends(get_db)):
    """
    Загружает JSON-файл с молекулами в БД.
    Ожидаемый формат файла:
    [{"structure": "CCO"}, {"structure": "CCC"}]
    """
    try:
        if file.content_type != "application/json":
            raise HTTPException(status_code=400,
                                detail="Только JSON-файлы поддерживаются.")

        file_content = await file.read()
        data = json.loads(file_content.decode("utf-8"))

        if not isinstance(data, list):
            raise HTTPException(status_code=400,
                                detail="Файл должен содержать список объектов.")

        for molecule in data:
            if "structure" not in molecule:
                raise HTTPException(
                    status_code=400,
                    detail="Каждый объект должен содержать поле 'structure'.")

            # Проверяем, что structure - валидный SMILES
            mol = Chem.MolFromSmiles(molecule["structure"])
            if mol is None:
                raise HTTPException(
                    status_code=400,
                    detail=f"Некорректный SMILES: {molecule['structure']}")

            db_molecule = Molecule(structure=molecule["structure"])
            db.add(db_molecule)

        db.commit()
        return {"message": "Молекулы успешно загружены.", "count": len(data)}

    except json.JSONDecodeError:
        raise HTTPException(status_code=400, detail="Некорректный формат JSON.")
