import pytest
from gem import matrix




def test_function():
    b=[[1.0, 0.0, 0.0, 0.0, 0.0],
 [0.0, 1.0, 0.0, 0.0, 0.0],
 [0.0, 0.0, 1.0, 0.0, 0.0],
 [0.0, 0.0, 0.0, 1.0, 0.0],
 [0.0, 0.0, 0.0, 0.0, 1.0]]
    a=matrix.identity(5)
    
    assert matrix.transpose(a)==b
def test_function1():
    b=[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    a=matrix.identity(3)
    
    assert matrix.transpose(a)==b
def test_function2():
    b=[[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
 [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
 [0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
 [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
 [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
 [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
 [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]]
    a=matrix.identity(7)
    
    assert matrix.transpose(a)==b





def test_function3():
  
    a=matrix.identity(3)
    
    assert matrix.inverse3(a)==a
def test_function5():
   
    a=matrix.identity(4)
    
    assert matrix.inverse4(a)==a

    