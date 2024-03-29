cd "$ROOT_DIR"/"$APP_DIR"
nohup redis-server > "$DATA_DIR"/redis-server.log 2>&1&
nohup celery -A "$APP_DIR" worker --loglevel=INFO --concurrency="${my_threads}" > "$DATA_DIR"/celery.log 2>&1&
python getSecretKey.py "$ROOT_DIR"/"$APP_DIR"
nohup python manage.py runserver 0.0.0.0:9000 > "$DATA_DIR"/django.log 2>&1&
cd "$DATA_DIR"
