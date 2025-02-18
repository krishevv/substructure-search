import json
from io import BytesIO
from fastapi.testclient import TestClient
from main import app

client = TestClient(app)


def test_add_molecule():
    response = client.post("/add", json={"id": 1, "smiles": "CCO"})
    assert response.status_code == 200
    assert response.json() == {
        "message": "Молекула добавлена успешно.", "id": 1}


def test_get_molecule():
    response = client.get("/molecule/1")
    assert response.status_code == 200
    assert response.json() == {"id": 1, "smiles": "CCO"}


def test_update_molecule():
    response = client.put("/molecule/1", json={"id": 1, "smiles": "CCC"})
    assert response.status_code == 200
    assert response.json() == {
        "message": "Молекула обновлена успешно.", "id": 1}


def test_delete_molecule():
    response = client.delete("/molecule/1")
    assert response.status_code == 200
    assert response.json() == {"message": "Молекула удалена успешно.", "id": 1}


def test_search_substructure():
    client.post("/add", json={"id": 2, "smiles": "CCO"})
    response = client.post("/search", json={"substructure": "CO"})
    assert response.status_code == 200
    assert response.json()["matches"] == [{"id": 2, "smiles": "CCO"}]


def test_add_duplicate_molecule():
    client.post("/add", json={"id": 3, "smiles": "CCC"})  # Добавляем молекулу
    # Пробуем добавить с тем же ID
    response = client.post("/add", json={"id": 3, "smiles": "CCO"})

    assert response.status_code == 400
    assert response.json()["detail"] == "Молекула с таким ID уже существует."


def test_get_nonexistent_molecule():
    response = client.get("/molecule/999")  # ID, которого нет

    assert response.status_code == 404
    assert response.json()["detail"] == "Молекула не найдена."


def test_upload_molecules():
    data = json.dumps([{"id": 10, "smiles": "CCN"},
                      {"id": 11, "smiles": "OCC"}])
    files = {
        "file": (
            "molecules.json",
            BytesIO(
                data.encode()),
            "application/json")}

    response = client.post("/upload", files=files)

    assert response.status_code == 200
    assert response.json() == {
        "message": "Молекулы успешно загружены.",
        "count": 2}
