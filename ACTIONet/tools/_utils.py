import string
import random

def rand_suffix(N):
    str_out = ''.join(random.choices(string.ascii_uppercase + string.ascii_lowercase+ string.digits, k=N))
    return str_out
