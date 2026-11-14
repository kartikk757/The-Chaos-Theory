"""Small process-local rate limiter for optional public archive clients."""
from __future__ import annotations

import threading
import time

_lock = threading.Lock()
_last_request: dict[str, float] = {}


def wait_for_service(service: str, minimum_interval_seconds: float = 0.5) -> None:
    """Ensure requests to one public service are spaced within this process."""
    interval = max(0.0, float(minimum_interval_seconds))
    with _lock:
        now = time.monotonic()
        previous = _last_request.get(service, 0.0)
        delay = max(0.0, interval - (now - previous))
        if delay:
            time.sleep(delay)
        _last_request[service] = time.monotonic()
