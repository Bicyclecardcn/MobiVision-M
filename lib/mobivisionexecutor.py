#!/usr/bin/env python
# ! -*- utf-8 -*-
"""
Analysis_Pipeline_Logging_System Command Executor Module: mobivisionexecutor.py
Author: jingxinxing
Date: 2025-05-06
Version: 1.0
"""
import os
import sys
import subprocess
import threading
from pathlib import Path
from datetime import datetime
from typing import List, Optional, Callable
from mobivisionlogging import MobiCommandLogSystem

class CommandExecutor:
    """支持双模式日志的命令执行器"""

    def __init__(
            self,
            log_system: MobiCommandLogSystem,
            console_output: bool = True, 
            #o_dir: str = None, 
            #dev_mode: bool = None, 
    ):
        self.log_system = log_system
        self.console_output = console_output
        self.logger = log_system.get_logger("command_executor")
        self.std_logs = log_system.std_logs

        # 初始化标准输出日志
        for log_file in self.std_logs.values():
            Path(log_file).touch(exist_ok=True)

    def execute(
            self,
            command: List[str],
            context: Optional[dict] = None,
            #callback: Optional[Callable] = None, 
            console_output: bool = True
    ) -> int:
        """执行命令并记录输出"""
        # full_cmd = " ".join(command)
        full_cmd = str(command)
        # self.logger.debug(
        #     f"[Executing Command] {full_cmd[0:20]} ...",
        #     command=full_cmd,
        #     context=context,
        #     env=self.log_system.env
        # )
        #self.logger.log("INFO", full_cmd)
        proc = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding="utf-8",
            errors="replace",
            shell=True
        )

        # 启动日志线程
        stdout_thread = threading.Thread(
            target=self._handle_stream,
            args=(proc.stdout, "stdout", context, console_output)
        )
        stderr_thread = threading.Thread(
            target=self._handle_stream,
            args=(proc.stderr, "stderr", context, console_output)
        )

        stdout_thread.start()
        stderr_thread.start()
        exit_code = -1
        try:
            exit_code = proc.wait()
        except subprocess.TimeoutExpired:
            proc.kill()
            self.logger.error(
                "[Command Timed Out]",
                command=full_cmd
            )
            raise
        finally:
            stdout_thread.join()
            stderr_thread.join()

        # 记录最终状态
        # log_method = self.logger.debug if exit_code == 0 else self.logger.error
        # log_method(
        #     f"[Command Completed] {full_cmd[0:20]} ...",
        #     command=full_cmd,
        #     exit_code=exit_code,
        #     context=context
        # )
        if exit_code != 0 or console_output:
            line = " cmd | " + full_cmd
            # 记录到文件
            with open(self.std_logs["stdout"], "a") as log_file:
                log_file.write(f"{datetime.now().isoformat()} | {line}\n")
        return exit_code

    def _handle_stream(self, stream, stream_type: str, context: dict, console_output: bool):
        """处理输出流"""
        with open(self.std_logs[stream_type], "a") as log_file:
            for line in iter(stream.readline, ""):
                line = line.strip()
                if not line:
                    continue
                line = " %s | " %(str(context))+ line
                # 记录到文件
                log_file.write(f"{datetime.now().isoformat()} | {line}\n")

                # 结构化日志记录
                log_context = {
                    "stream": stream_type,
                    "content": line,
                    "context": context
                }

                #if console_output:
                #    if self.log_system.env == "Production":
                #        level = "INFO" if stream_type == "stdout" else "ERROR"
                #        self.logger.log(level, line)
                #    else:
                #        self.logger.debug("[Command Output: ]", line)

                # 控制台输出
                #if console_output:
                #    prefix = "[STDOUT]" if stream_type == "stdout" else "[STDERR]"
                #    print(f"{prefix} {line}")