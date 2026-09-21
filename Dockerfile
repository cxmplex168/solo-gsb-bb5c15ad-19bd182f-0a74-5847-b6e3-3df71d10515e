FROM python:3.12-slim

ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1

WORKDIR /app

# 依赖安装（构建镜像时完成）
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# 应用代码
COPY app ./app
COPY tests ./tests

# SQLite 数据目录（可挂载卷持久化）
RUN mkdir -p /data
ENV APP_DB_PATH=/data/digest.db \
    APP_HOST=0.0.0.0 \
    APP_PORT=8000

EXPOSE 8000

# 启动即初始化数据库（应用 lifespan 中建表），随后对外服务
CMD ["sh", "-c", "uvicorn app.main:app --host ${APP_HOST} --port ${APP_PORT}"]
