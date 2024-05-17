FROM python:3
WORKDIR /app
COPY . .
RUN pip install --no-cache-dir -r requirements.txt
RUN tar -xf references/chr6.tar.xz
CMD ["python", "./main.py"]

