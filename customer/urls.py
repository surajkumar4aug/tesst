# from django.urls import path
# from .views import custermerupdate,User
# from django.contrib import admin
# from rest_framework.urlpatterns import format_suffix_patterns

# urlpatterns = [
   
    
#     path('userdetails/<int:pk>', custermerupdate.as_view(), name='userdetail'),
#     path('userdetails/', User.as_view(), name='user')
# ]
# from django.urls import path
# from .views import AuthorList, BookList, CategoryList

# urlpatterns = [
#     path('authors/', AuthorList.as_view(), name='author-list'),
#     path('books/', BookList.as_view(), name='book-list'),
#     path('categories/', CategoryList.as_view(), name='category-list'),
# ]

from django.urls import path
from .views import AuthorAPIView, BookAPIView, PublisherAPIView
from .views import MusicView

urlpatterns = [
    path('authors/', AuthorAPIView.as_view(), name='author-list'),
    path('authors/<int:pk>/', AuthorAPIView.as_view(), name='author-detail'),
    path('books/', BookAPIView.as_view(), name='book-list'),
    path('books/<int:pk>/', BookAPIView.as_view(), name='book-detail'),
    path('publishers/', PublisherAPIView.as_view(), name='publisher-list'),
    path('publishers/<int:pk>/', PublisherAPIView.as_view(), name='publisher-detail'),
    path('music/', MusicView.as_view(), name='music-list'),
]



