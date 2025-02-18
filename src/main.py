from fastapi import UploadFile, File, HTTPException
from fastapi import FastAPI, HTTPException, File
from pydantic import BaseModel
from rdkit import Chem
from os import getenv

import json


# Инициализация приложения FastAPI
app = FastAPI()

# Глобальное хранилище молекул
molecule_store = {}

# Модели данных


class Molecule(BaseModel):
    id: int
    smiles: str


class SubstructureQuery(BaseModel):
    substructure: str

# Эндпоинты


@app.get("/")
def get_server():
    return {"server_id": getenv("SERVER_ID", "unknown")}


@app.post("/add")
def add_molecule(molecule: Molecule):
    if molecule.id in molecule_store:
        raise HTTPException(status_code=400,
                            detail="Молекула с таким ID уже существует.")
    molecule_store[molecule.id] = molecule.smiles
    return {"message": "Молекула добавлена успешно.", "id": molecule.id}


@app.get("/molecule/{id}")
def get_molecule(id: int):
    if id not in molecule_store:
        raise HTTPException(status_code=404, detail="Молекула не найдена.")
    return {"id": id, "smiles": molecule_store[id]}


@app.put("/molecule/{id}")
def update_molecule(id: int, molecule: Molecule):
    if id not in molecule_store:
        raise HTTPException(status_code=404, detail="Молекула не найдена.")
    molecule_store[id] = molecule.smiles
    return {"message": "Молекула обновлена успешно.", "id": id}


@app.delete("/molecule/{id}")
def delete_molecule(id: int):
    if id not in molecule_store:
        raise HTTPException(status_code=404, detail="Молекула не найдена.")
    del molecule_store[id]
    return {"message": "Молекула удалена успешно.", "id": id}


@app.get("/list")
def list_molecules():
    return {"molecules": molecule_store}


@app.post("/search")
def search_substructure(query: SubstructureQuery):
    substructure = Chem.MolFromSmiles(query.substructure)
    if not substructure:
        raise HTTPException(
            status_code=400,
            detail="Некорректная SMILES строка для подструктуры.")

    matching_molecules = []
    for id, smiles in molecule_store.items():
        mol = Chem.MolFromSmiles(smiles)
        if mol and mol.HasSubstructMatch(substructure):
            matching_molecules.append({"id": id, "smiles": smiles})

    return {"matches": matching_molecules}


@app.post("/upload")
async def upload_file(file: UploadFile = File(...)):
    """
    Эндпоинт для загрузки файла с молекулами.
    Ожидает JSON-файл с массивом объектов вида:
    [{"id": 1, "smiles": "CCO"}, {"id": 2, "smiles": "CCC"}]
    """
    try:
        # Проверяем, что это JSON-файл
        if file.content_type != "application/json":
            raise HTTPException(status_code=400,
                                detail="Только JSON-файлы поддерживаются.")

        # Читаем содержимое файла
        file_content = await file.read()
        data = json.loads(file_content.decode("utf-8"))

        # Проверяем содержимое файла
        if not isinstance(data, list):
            raise HTTPException(
                status_code=400,
                detail="Файл должен содержать список объектов.")

        # Добавляем молекулы в наш контейнер
        for molecule in data:
            if "id" not in molecule or "smiles" not in molecule:
                raise HTTPException(
                    status_code=400,
                    detail="Каждый объект должен содержать поля 'id' и 'smiles'.")
            if molecule["id"] in molecule_store:
                raise HTTPException(
                    status_code=400, detail=f"Молекула с ID \
                          {molecule['id']} уже существует.")
            molecule_store[molecule["id"]] = molecule["smiles"]

        return {"message": "Молекулы успешно загружены.", "count": len(data)}

    except json.JSONDecodeError:
        raise HTTPException(
            status_code=400,
            detail="Некорректный формат JSON.")
