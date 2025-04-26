**Substructure Search**  — это веб-приложение для поиска подструктур в химических соединениях. Оно использует Python и Nginx, развертывается с помощью Docker Compose.​

##  Структура проекта 

 
- `src/` — исходный код приложения.
 
- `nginx/` — конфигурационные файлы Nginx.
 
- `docker-compose.yml` — файл для развертывания проекта с использованием Docker Compose.​


## Быстрый старт 

 
2. Убедитесь, что у вас установлены [Docker](https://www.docker.com/)  и [Docker Compose](https://docs.docker.com/compose/) .
 
4. Клонируйте репозиторий:​


```bash
git clone https://github.com/krishevv/substructure-search.git
cd substructure-search
```


3. Запустите приложение:​


```bash
docker-compose up --build
```


4. Откройте веб-браузер и перейдите по адресу `http://localhost` для доступа к приложению.​


 
