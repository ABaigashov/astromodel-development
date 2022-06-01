import logging.config
from logging_config import LOGGING_CONFIG


class Logger:
    def __init__(self):
        logging.config.dictConfig(LOGGING_CONFIG)
        self.logger = logging.getLogger('my_logger')
