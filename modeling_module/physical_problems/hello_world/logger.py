
class CreationException(Exception):
    pass


class DataError(CreationException):
    def __init__(self):  # custom property
        super().__init__(
            f"Data"
        )

