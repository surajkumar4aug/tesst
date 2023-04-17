# from xml.parsers.expat import model
# from django.db import models
# from django.forms import CharField, EmailField

# class userDetail(models.Model):
#     name=models.CharField(max_length=10,blank=True, null=True)
#     email_id=models.EmailField(blank=True, null=True)
#     phone_no = models.M(max_length = 10,blank=True, null=True)
#     address = models.CharField(max_length = 50,blank=True, null=True)
#     class meta:
#         db_table = 'userDetails'
#     def __str__(self):
# 		    return self.name
# class userDetail(models.Model):
#     name=models.CharField(max_length=10,blank=True, null=True)
#     email_id=models.EmailField(blank=True, null=True)
#     phone_no = models.CharField(max_length = 10,blank=True, null=True)
#     address = models.CharField(max_length = 50,blank=True, null=True)
#     class meta:
#         db_table = 'userDetails'
#     def __str__(self):
# 		    return self.name    
# from django.db import models

# class Author(models.Model):
#     name = models.CharField(max_length=255)

#     def __str__(self):
#         return self.name

# class Book(models.Model):
#     title = models.CharField(max_length=255)
#     author = models.ForeignKey(Author, on_delete=models.CASCADE)

#     def __str__(self):
#         return self.title

# class Category(models.Model):
#     name = models.CharField(max_length=255)
#     books = models.ManyToManyField(Book)

#     def __str__(self):
#         return self.name    




from django.db import models

class Publisher(models.Model):
    name = models.CharField(max_length=100)
    website = models.URLField()

    def __str__(self):
        return self.name

class Author(models.Model):
    name = models.CharField(max_length=100)
    email = models.EmailField(null=True)
    #books = models.ManyToManyField('Book', related_name='authors')

    def __str__(self):
        return self.name

class Book(models.Model):
    title = models.CharField(max_length=100)
    description = models.TextField(null=True)
    publisher = models.ForeignKey(Publisher,null=True, on_delete=models.CASCADE)
    authors = models.ManyToManyField(Author)

    def __str__(self):
        return self.title
from django.db import models

class Artist(models.Model):
    name = models.CharField(max_length=100)
    # ...

class Album(models.Model):
    title = models.CharField(max_length=100)
    artist = models.ForeignKey(Artist, on_delete=models.CASCADE)
    # ...

class Song(models.Model):
    title = models.CharField(max_length=100)
    artist = models.ForeignKey(Artist, on_delete=models.CASCADE)
    album = models.ForeignKey(Album, on_delete=models.CASCADE)
    # ...
from django.db import models
from django.contrib.auth.models import User

class Chat(models.Model):
    sender = models.ForeignKey(User, related_name='sender', on_delete=models.CASCADE)
    recipient = models.ForeignKey(User, related_name='recipient', on_delete=models.CASCADE)
    message = models.TextField()
    timestamp = models.DateTimeField(auto_now_add=True)