#threads=8
threads=16
#threads=4
cd "$ROOT_DIR"/"$APP_DIR"
nohup redis-server > "$DATA_DIR"/redis-server.log 2>&1&
nohup celery -A "$APP_DIR" worker --loglevel=INFO --concurrency="$threads" > "$DATA_DIR"/celery.log 2>&1&
nohup python manage.py runserver 0.0.0.0:9000 > "$DATA_DIR"/django.log 2>&1&
cd "$DATA_DIR"
