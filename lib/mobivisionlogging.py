#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""
Analysis Pipeline Logging System: mobivisionlogging.py, mobivisionexecutor.py, main.py
Author: jingxinxing
Description:
    (1)This module provides a structured logging configuration for the bioinformatics analysis pipeline.
    (2)It includes a structured logging configuration that can be used across the project.
    (3)It's support for both production and development environments.
    (4)It's has mutliple logging system, one for Normal User, one for Develop  User.
    (5)It's support asyncio logging system.

Version: 1.0
Date: 2025-03-27
Update: 2025-04-23
"""

## Import Modules
import os
import sys
import json
import rtoml
import logging
import asyncio
import structlog
import threading
import subprocess
import logging.config
from pathlib import Path
from datetime import datetime
from loguru import logger as user_logger
from typing import Dict, Any, Optional, Literal

## Type Hints
# -------------------------------#
LogConfig = Dict[str, Any]

# Config_File_Path = "/PATH/TO/LOGCONFIG"

LogEnvType = Literal["Development", "Production"]

LogLevelType = Literal["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]
# ------------------------------#

def filter_exclude_star(record):
    # 定义关键词列表和截取长度
    STAR_KEYWORDS = [
        "STAR --runThreadN",
        "STAR version",
        "started STAR run",
        "loading genome",
        "mapping",
        "Solo counting",
        "sorting BAM",
        "finished successfully"
    ]
    message = record["message"]
    return not any(keyword in message for keyword in STAR_KEYWORDS)

def filter_exclude_command(record):
    # 过滤关键词
    COMMAND_KEYWORDS = [
        "Command",
        "STAR --runThreadN",
        "STAR version",
        "started STAR run",
        "loading genome",
        "mapping",
        "Solo counting",
        "sorting BAM",
        "finished successfully"
    ]
    message = record["message"]
    return not any(keyword in message for keyword in COMMAND_KEYWORDS)

## Class: MobiLoggingSystem
class MobiLoggingSystem:
    """MobiVision Logging Module"""
    def __init__(self, o_dir: str = None, 
                 dev_mode: bool = None, 
                 config_path: str = None,
                 log_message: str = None,
                 print_level: LogLevelType = "INFO",
                 logfile_level: LogLevelType = "DEBUG",
                 log_level: LogLevelType = "DEBUG"
                 ): # command: str = None
        if not o_dir is None:
            self.working_path = o_dir
        else:
            if "-o" in sys.argv[1:] or "--outDir" in sys.argv[1:] or "--resultDir" in sys.argv[1:]:
                self.working_path = Path(sys.argv[sys.argv.index("-o") + 1])
            elif "--outDir" in sys.argv[1:]:
                self.working_path = Path(sys.argv[sys.argv.index("--outDir") + 1])
            elif "--resultDir" in sys.argv[1:]:
                self.working_path = Path(sys.argv[sys.argv.index("--resultDir") + 1])
            else:
                self.working_path = Path.cwd()
        self.working_path = str(self.working_path)
        ## Developer Signals
        if not dev_mode is None:
            if dev_mode:
                self.env = "Development"
            else:
                self.env = "Production"
            self.dev_mode = dev_mode
        else:
            if (Path.cwd() / "Development").exists():
                self.env = "Development"
                self.dev_mode = True
            else:
                self.env = "Production"
                self.dev_mode = False
        self.config_path = config_path
        # self.env = env
        self.log_message = log_message
        self.print_level = print_level
        self.logfile_level = logfile_level
        # self.module_name = module_name
        self.log_level = log_level
        # self.project_root_path = project_root_path
        # self.command = command

    ## Function: _get_mobivision_version
    @classmethod
    def _get_mobivision_version(cls):
        """
        Get the version of the MobiVision software.
        Args:
            mobivision_help_file (str): The path to the MobiVision help file.
        Returns:
            Software_Version, str: The version of the MobiVision software.
        """
        mobivision_help_dir = Path(__file__).resolve().parent.resolve().parent / "config"
        mobivision_help_file = mobivision_help_dir / "mobivision_help.toml"

        with open(mobivision_help_file, "r") as f:
            config = rtoml.load(f)
            for section in config:
                content = config[section]
                if "mobivision_version" in content:
                    cls.Software_Version = content["mobivision_version"]
                    break
        return cls.Software_Version

    ## Function: _mobilogrecorder(_mobivision_logging_recorder -> _mobilogrecorder)
    # def _mobilogrecorder(self, project_root_path, log_message, module_name, log_level, env, config_path: str = None, print_level = "INFO", logfile_level = "DEBUG"):
    def _mobilogrecorder(self, log_message, log_level, config_path: str = None,
                         print_level="INFO", logfile_level="DEBUG", show = True):
        """
        根据运行环境(env)控制控制台日志输出，同时确保所有日志写入文件。
        Args:
            project_root_path: 项目根目录路径
            log_message: 日志消息
            module_name: 模块名称
            log_level: 日志级别
            config_path: 日志配置文件路径
            print_level: 控制台日志级别
            logfile_level: 日志文件日志级别
        """

        ## User Logger ##
        log_name = "MobiVision-M"

        # 根据参数进行执行对应的日志功能语句：env、log_level、log_message
        if self.env not in ("Development", "Production"):
            raise ValueError(f"Invalid env: {self.env}, must be 'Development' or 'Production'")
        else:
            # print(f"MobiVision Running in {env} Mode.")
            if self.env == "Development":
                #
                ## Developer Logger ##
                if not config_path:
                    module_path = Path(__file__).resolve().parent.resolve().parent
                    config_path = module_path / "config" / "logging_config.json"
                # Config_File_Path = config_path
                # 加载基础配置
                with open(config_path, "r", encoding="utf-8") as f:  # Here
                    config: LogConfig = json.load(f)

                # 动态生成当前工作目录中的日志路径
                log_dir = self.working_path / "logs" # project_root_path -> working_path
                log_dir.mkdir(parents=True, exist_ok=True)

                # 更新文件处理器路径
                config["handlers"]["file"]["filename"] = str(log_dir / f".{log_name}.log")

                # 应用配置
                logging.config.dictConfig(config)

                # 结构化日志配置
                structlog.configure(
                    processors=[
                        structlog.stdlib.filter_by_level,
                        structlog.stdlib.add_logger_name,
                        structlog.stdlib.add_log_level,
                        structlog.processors.TimeStamper(fmt="iso"),
                        structlog.processors.JSONRenderer()
                    ],
                    context_class=dict,
                    logger_factory=structlog.stdlib.LoggerFactory(),
                    wrapper_class=structlog.stdlib.BoundLogger,
                    cache_logger_on_first_use=True,  # True, False
                )

                # 获取模块级Logger
                developer_logger = structlog.get_logger(__name__)
                developer_logger.propagate = False
                #

                if log_level not in ("DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"):
                    raise ValueError(
                        f"Invalid log_level: {log_level}, must be 'DEBUG', 'INFO', 'WARNING', 'ERROR', or 'CRITICAL'")
                else:
                    # print(f"MobiVision _mobivision_logging_recorder() Logging Level: {log_level}.")
                    if log_level == "DEBUG":
                        developer_logger.debug(log_message)
                        # user_logger.debug(log_message)
                    elif log_level == "INFO":
                        developer_logger.info(log_message)
                        # user_logger.info(log_message)
                    elif log_level == "WARNING":
                        developer_logger.warning(log_message)
                        # user_logger.warning(log_message)
                    elif log_level == "ERROR":
                        developer_logger.error(log_message)
                        # user_logger.error(log_message)
                    elif log_level == "CRITICAL":
                        developer_logger.critical(log_message)
                        # user_logger.critical(log_message)
                    else:
                        raise ValueError(
                            f"Invalid log_level: {log_level}, must be 'DEBUG', 'INFO', 'WARNING', 'ERROR', or 'CRITICAL'")
            elif self.env == "Production":
                # 用户日志记录器 (User Logger)
                user_logger.remove()
                # user_logger.add(sys.stdout, level=print_level, format="{time:YYYY-MM-DD HH:mm:ss} | <green>{level}</green> | <blue>{message}</blue>",
                #                 colorize=True, enqueue=True, backtrace=True, diagnose=False) # filter=lambda record: env == "Production",  enqueue=True, backtrace=True, diagnose=True)
                # user_logger.add(sys.stdout, level=print_level, format=custom_formatter, colorize=True, enqueue=True, backtrace=True, diagnose=False)
                if show:
                    user_logger.add(sys.stdout, level=print_level,
                                    format="{time:YYYY-MM-DD HH:mm:ss} | <green>{level: ^10}</green> | <blue>{message}</blue>",
                                    colorize=True, enqueue=True, backtrace=True, diagnose=False, filter=filter_exclude_star)
                # user_logger.add(sys.stderr, level="ERROR",
                #                 format="{time:YYYY-MM-DD HH:mm:ss} | <red>{level}</red> | <blue>{message}</blue>",
                #                 colorize=True,
                #                 enqueue=True, backtrace=True, diagnose=True)
                user_logger.add(sink=os.path.join(self.working_path, "logs", f"{log_name}.log"),
                                format="{time:YYYY-MM-DD HH:mm:ss} | {level: ^10} | {message}", level=logfile_level,
                                rotation="100 MB", enqueue=True, backtrace=True, diagnose=True,
                                filter=filter_exclude_command)
                #
                if log_level not in ("DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"):
                    raise ValueError(f"Invalid log_level: {log_level}, must be 'DEBUG', 'INFO', 'WARNING', 'ERROR', or 'CRITICAL'")
                else:
                    # print(f"MobiVision _mobivision_logging_recorder() Logging Level: {log_level}.")
                    if log_level == "DEBUG":
                        # developer_logger.debug(log_message)
                        user_logger.debug(log_message)
                    elif log_level == "INFO":
                        # developer_logger.info(log_message)
                        user_logger.info(log_message)
                    elif log_level == "WARNING":
                        # developer_logger.warning(log_message)
                        user_logger.warning(log_message)
                    elif log_level == "ERROR":
                        # developer_logger.error(log_message)
                        user_logger.error(log_message)
                    elif log_level == "CRITICAL":
                        # developer_logger.critical(log_message)
                        user_logger.critical(log_message)
                    else:
                        raise ValueError(f"Invalid log_level: {log_level}, Must be 'DEBUG', 'INFO', 'WARNING', 'ERROR', or 'CRITICAL'")
            else:
                raise ValueError(f"Invalid env: {self.env}, Must be 'Development' or 'Production'")
        return ""

## Class: MobiCommandLogSystem (MobiCommandLogSystem -> MobiCommandLogSystem)
class MobiCommandLogSystem:
    """双模式日志系统核心类"""

    def __init__(
            self,
            o_dir: str = None, 
            dev_mode: bool = None, 
            # working_path: Path,
            # env: Literal["Production", "Development"] = "Production",
            log_level: str = "DEBUG" # INFO -> DEBUG
    ):
        if not o_dir is None:
            self.working_path = Path(o_dir)
        else:
            if "-o" in sys.argv[1:] or "--outDir" in sys.argv[1:] or "--resultDir" in sys.argv[1:]:
                self.working_path = Path(sys.argv[sys.argv.index("-o") + 1])
            elif "--outDir" in sys.argv[1:]:
                self.working_path = Path(sys.argv[sys.argv.index("--outDir") + 1])
            elif "--resultDir" in sys.argv[1:]:
                self.working_path = Path(sys.argv[sys.argv.index("--resultDir") + 1])
            else:
                self.working_path = Path.cwd()

        ## Developer Signals
        if not dev_mode is None:
            if dev_mode:
                self.env = "Development"
            else:
                self.env = "Production"
        else:
            if (Path.cwd() / "Development").exists():
                self.env = "Development"
            else:
                self.env = "Production"
        self.log_level = log_level
        self.log_dir = self.working_path / "logs"
        self._configure()

    def _configure(self):
        """配置日志系统"""
        self.log_dir.mkdir(parents=True, exist_ok=True)

        if self.env == "Production":
            self._setup_production()
        else:
            self._setup_development()

    def _setup_production(self):
        """生产环境配置"""
        # 移除默认配置
        user_logger.remove()

        # 控制台配置：外部程序运行时关闭控制台输出
        user_logger.add(
            sink=sys.stdout,
            level=self.log_level,
            format="<green>{time:YYYY-MM-DD HH:mm:ss}</green> | <level>{level: <8}</level> | {message}",
            colorize=True,
            backtrace=True,
            diagnose=False,
            filter=filter_exclude_star
        )
        # 主日志文件
        user_logger.add(
            sink=self.log_dir / "MobiVision-M.log",
            level=self.log_level,
            rotation="100 MB",
            format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}",
            enqueue=True,
            backtrace=True,
            filter=filter_exclude_command
        )

    def _setup_development(self):
        """开发环境配置"""
        # 结构化日志配置
        structlog.configure(
            processors=[
                structlog.contextvars.merge_contextvars,
                structlog.stdlib.filter_by_level,
                structlog.stdlib.add_logger_name,
                structlog.stdlib.add_log_level,
                structlog.processors.TimeStamper(fmt="iso"),
                structlog.processors.StackInfoRenderer(),
                structlog.processors.format_exc_info,
                structlog.processors.JSONRenderer()
            ],
            wrapper_class=structlog.stdlib.BoundLogger,
            context_class=dict,
            logger_factory=structlog.stdlib.LoggerFactory(),
            cache_logger_on_first_use=False,
        )

        # 标准logging配置
        logging.basicConfig(
            format="%(asctime)s | %(levelname)s | %(module)s | %(funcName)s | %(lineno)d | %(message)s",
            level=logging.DEBUG,
            handlers=[
                logging.FileHandler(self.log_dir / ".MobiVision-M.log"),
                logging.StreamHandler()
            ]
        )

    def get_logger(self, name: str = None) -> Any:
        """获取对应环境的日志记录器"""
        if self.env == "Production":
            return user_logger
        return structlog.get_logger(name or __name__)

    @property
    def std_logs(self) -> Dict[str, Path]:
        """获取标准输出日志路径"""
        return {
            "stdout": os.path.join(self.log_dir, "stdout.log"),
            "stderr": os.path.join(self.log_dir, "stderr.log")
        }