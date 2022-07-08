class MyError(Exception):
    def __init__(self, text):
        self.txt = text


class Logger:
    def __init__(self, general_array, objects_array):
        self.print_general(general_array)

    @classmethod
    def print_general(cls, general_array):
        print(general_array)
        for pair in general_array:
            field = pair['field']
            error = pair['error']
            print(f'Error { error } in { field }') #
            

class Tester:
    def __init__(self, config):
        self.general_result = self.test_general(config._CFG_RAW_DATA['general'])
        self.objects_result = self.test_objects(config._CFG_RAW_DATA['objects'])

        self.log = Logger(self.general_result, self.objects_result)

    @classmethod
    def test_general(cls, general):
        result = []
        for field in general:
            if type(general[field]) == "<class 'dict'>": #
                value = general[field]['value']
            value = general[field]

            if cls.test_empty(value):
                result += [{'field': field, 'error': 'empty value'}]
            
            
        return result

    @classmethod
    def test_objects(cls, objects):
        pass

        

    @classmethod
    def get_objectValue(cls, object):
        pass
        # return object['value']

    @classmethod
    def test_type(cls, value):
        pass

    @classmethod
    def test_empty(cls, value):
        return value == ""