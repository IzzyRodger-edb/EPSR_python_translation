"""The shell on a real terminal (pexpect): prompt, built-ins, Ctrl+C, shutdown, no residual
processes (PLAN.md section 4)."""

import os
import signal
import sys
import time

import pexpect
import psutil
import pytest

PROMPT = "EPSRshell> "


def spawn(folder):
    env = dict(os.environ, PYTHONUNBUFFERED="1", TERM="dumb")
    child = pexpect.spawn(sys.executable, ["-m", "epsr26"], cwd=str(folder), env=env, encoding="utf-8", timeout=20)
    child.expect("Welcome to EPSR version 26")
    child.expect(PROMPT)
    return child


def cmd(child, line):
    child.sendline(line)
    child.expect(PROMPT)
    return child.before


def no_processes_like(pattern):
    for p in psutil.process_iter(["cmdline"]):
        try:
            if pattern in " ".join(p.info["cmdline"] or []):
                return False
        except psutil.NoSuchProcess:
            pass
    return True


def wait_exit(child, timeout=10):
    child.expect(pexpect.EOF, timeout=timeout)
    child.close()
    return child.exitstatus, child.signalstatus


def test_builtins(tmp_path):
    (tmp_path / "sub").mkdir()
    (tmp_path / "sub" / "afile.txt").write_text("x\n")
    c = spawn(tmp_path)
    out = cmd(c, "pwd")
    assert f"Current folder for input files: {tmp_path}/" in out
    out = cmd(c, "cd sub")
    assert f"Current folder for input files: {tmp_path}/sub/" in out
    out = cmd(c, "ls")
    assert "afile.txt" in out
    out = cmd(c, "ls -l afile.txt")
    assert "afile.txt" in out and "-rw" in out
    out = cmd(c, "cd ..")
    assert f"Current folder for input files: {tmp_path}/" in out
    out = cmd(c, "cd -")
    assert f"{tmp_path}/sub/" in out
    out = cmd(c, "cd /no/such/folder")
    assert "No such folder!" in out
    out = cmd(c, "!echo hi-from-sh")
    assert "hi-from-sh" in out
    out = cmd(c, "md newdir")
    assert (tmp_path / "sub" / "newdir").is_dir()
    out = cmd(c, "del afile.txt")
    assert not (tmp_path / "sub" / "afile.txt").exists()
    c.sendline("bye")
    status, _ = wait_exit(c)
    assert status == 0
    assert not (tmp_path / ".epsrshell.lock").exists()
    assert (tmp_path / "runflag.txt").read_text().strip() == "interactive"


def test_ctrl_c_at_prompt_and_ctrl_d(tmp_path):
    c = spawn(tmp_path)
    c.send("half typed line")
    c.sendintr()
    c.expect(PROMPT)
    out = cmd(c, "pwd")  # the shell is still alive and the discarded line was not run
    assert "Current folder" in out
    c.sendeof()
    status, _ = wait_exit(c)
    assert status == 0


def test_no_residual_processes_after_bye(tmp_path):
    tag = f"sleep {300 + os.getpid() % 100}.5"
    c = spawn(tmp_path)
    cmd(c, f"system {tag} &")
    out = cmd(c, "jobs")
    assert tag in out
    assert not no_processes_like(tag)
    c.sendline("bye")
    wait_exit(c)
    time.sleep(0.5)
    assert no_processes_like(tag)


def test_sighup_shuts_down_cleanly(tmp_path):
    tag = f"sleep {400 + os.getpid() % 100}.5"
    c = spawn(tmp_path)
    cmd(c, f"system {tag} &")
    os.kill(c.pid, signal.SIGHUP)
    c.expect(pexpect.EOF, timeout=15)
    c.close()
    time.sleep(0.5)
    assert no_processes_like(tag)
    assert not (tmp_path / ".epsrshell.lock").exists()


def test_lock_warning(tmp_path):
    a = spawn(tmp_path)
    b = pexpect.spawn(sys.executable, ["-m", "epsr26"], cwd=str(tmp_path), encoding="utf-8", timeout=20,
                      env=dict(os.environ, PYTHONUNBUFFERED="1", TERM="dumb"))
    b.expect("another EPSRshell is running in this folder: PID")
    b.expect(PROMPT)
    b.sendline("bye")
    wait_exit(b)
    a.sendline("bye")
    wait_exit(a)


@pytest.mark.parametrize("script", ["# comment\npwd\nes\n"])
def test_script_mode(tmp_path, script):
    (tmp_path / "s.txt").write_text(script)
    c = spawn(tmp_path)
    c.sendline("ss s.txt")
    c.expect("Running script file: s.txt")
    c.expect("Current folder for input files")
    c.expect(PROMPT)
    c.sendline("bye")
    wait_exit(c)
    assert (tmp_path / "runflag.txt").read_text().strip() == "interactive"
