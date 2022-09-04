import os
from dotenv import dotenv_values
import logging
import argparse

# TODO: Make training/test samples ratio customizable

class EnvironmentVariableNotFoundError(Exception):
    def __init__(self, var, allowed_vars):
        message = f'{var}.\nAllowed variables in the environment are: {allowed_vars}'
        super().__init__(message)


class InvalidRepresentationTypeError(Exception):
    def __init__(self, repr_type):
        message = f'{repr_type}.\nRepresentation type can only be "binary" or "coverage".'
        super().__init__(message)


class Environment:
    def __init__(self):
        self.variables = {
            'ROOT_DIR': os.getcwd(),
            'REPR_TYPE': 'coverage',
            'ENV_FILE': '.env',
            'ENV_INITIALIZED': True,
            'NUM_OF_REGIONS': 30
        }
        self.allowed_vars = [
            k for k in self.variables.keys() if k not in ['ENV_INITIALIZED']
        ]

    @classmethod
    def __init_from_env__(cls):
        env = cls()
        env.variables = dotenv_values('.env')
        env.allowed_vars = ['ROOT_DIR', 'REPR_TYPE', 
                            'ENV_FILE', 'NUM_OF_REGIONS']
        return env
    
    def update_env(self):
        path_to_env_file = os.path.join(
            self.variables['ROOT_DIR'], 
            self.variables['ENV_FILE']
        )
        with open(path_to_env_file, 'w') as env_file:
            for var, val in self.variables.items():
                env_file.write(f'{var} = {val}\n')
    
    def uninitialize_env(self):
        self.variables = {
            'ROOT_DIR': os.getcwd(),
            'ENV_FILE': '.env',
            'ENV_INITIALIZED': ''
        }
        self.update_env()

    def change_var(self, var, val):
        if var not in self.allowed_vars:
            raise EnvironmentVariableNotFoundError(var, self.allowed_vars)
        elif var == 'REPR_TYPE' and val not in ['binary', 'coverage']:
            raise InvalidRepresentationTypeError(val)
        elif var == 'ENV_FILE' and not os.path.exists(
            os.path.join(self.variables['ROOT_DIR'], val)
        ):
            raise FileNotFoundError(os.path.join(self.variables['ROOT_DIR'], val))
        
        self.variables[var] = val
        self.update_env()


if __name__ == '__main__':
    logging.getLogger().setLevel(logging.INFO)
    if not dotenv_values('.env')['ENV_INITIALIZED']:
        env = Environment()
        env.update_env()
        logging.info('No initialized environment was found. Environment initialized.')
    else:
        env = Environment.__init_from_env__()
        logging.info('A previously initialized environment was found. Variables loaded.')
    parser = argparse.ArgumentParser(
        description='Parser for setup CLI'
    )
    parser.add_argument('--uninitialize', action='store_true')

    args = parser.parse_args()
    if args.uninitialize:
        logging.info('Uninitializing environment. ROOT_DIR variable will be kept.')
        env.uninitialize_env()
