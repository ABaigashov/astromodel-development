LOGGING_CONFIG = {
    'version': 1,
    'disable_existing_loggers': False,

    'formatters': {
        'default_formatter': {
            'format': '[%(levelname)s:%(asctime)s] %(message)s'
        },
    },

    'handlers': {
        'file_handler': {
            'class': 'logging.FileHandler',
            'filename': 'logconfig.log',
            'formatter': 'default_formatter',
        },
    },

    'loggers': {
        'my_logger': {
            'handlers': ['file_handler'],
            'level': 'DEBUG',
            'propagate': True
        }
    }
}