#!/usr/bin/env python3

import fcntl
import os
import pty
import struct
import subprocess
import sys
import termios
import time


def wait_for_quiet(process, master, transcript, deadline):
    last_data = None
    while process.poll() is None:
        now = time.monotonic()
        if now >= deadline:
            return False
        try:
            data = os.read(master, 65536)
            if data:
                transcript.extend(data)
                last_data = now
        except BlockingIOError:
            if last_data is not None and now - last_data >= 0.02:
                return True
            time.sleep(0.0005)
        except OSError:
            pass
    return False


def main():
    if len(sys.argv) < 2:
        return 2
    master, slave = pty.openpty()
    fcntl.ioctl(slave, termios.TIOCSWINSZ, struct.pack("HHHH", 24, 80, 0, 0))
    process = subprocess.Popen(sys.argv[1:], stdin=slave, stdout=slave, stderr=slave)
    os.close(slave)
    os.set_blocking(master, False)
    transcript = bytearray()
    deadline = time.monotonic() + 30
    if not wait_for_quiet(process, master, transcript, deadline):
        process.kill()
        process.wait()
        os.close(master)
        return 124
    started = time.monotonic_ns()
    os.write(master, b"l")
    if not wait_for_quiet(process, master, transcript, deadline):
        process.kill()
        process.wait()
        os.close(master)
        return 124
    elapsed = (time.monotonic_ns() - started) / 1_000_000_000
    os.write(master, b"q")
    while process.poll() is None:
        if time.monotonic() >= deadline:
            process.kill()
            process.wait()
            os.close(master)
            return 124
        try:
            transcript.extend(os.read(master, 65536))
        except BlockingIOError:
            time.sleep(0.001)
        except OSError:
            pass
    os.close(master)
    if process.returncode != 0:
        sys.stderr.buffer.write(transcript)
        return process.returncode
    print(f"{elapsed:.9f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
