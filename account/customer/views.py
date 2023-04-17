# from rest_framework.views import APIView
# from rest_framework.response import Response
# from .serializer import customerserial
# from .models import userDetail


# # Create your views here.
# class User(APIView):
#     def get(self, request):
#         userdata=userDetail.objects.all()
#         userSerializer=customerserial(userdata,many=True)
#         return Response(userSerializer.data)
#     def post(self, request):
#         serializer = customerserial(data=request.data)
#         if serializer.is_valid():
#             serializer.save()
#             return Response(serializer.data, status=201)
#         return Response(serializer.errors, status=400)
# class custermerupdate(APIView)   :
#     def get(self, request,pk):
#         data=userDetail.objects.get(pk=pk)
#         userSerializer=customerserial(data)
#         return Response(userSerializer.data)
    
        
#     def put(self,request,pk):
#         data=userDetail.objects.get(pk=pk)
#         serializer = customerserial(data,data=request.data)
#         if serializer.is_valid():
#             serializer.save()
#             return Response(serializer.data, status=201)
#         return Response(serializer.errors, status=400)     


# from rest_framework import generics
# from .models import Author, Book, Category
# from .serializer import AuthorSerializer, BookSerializer, CategorySerializer

# class AuthorList(generics.ListCreateAPIView):
#     queryset = Author.objects.all()
#     serializer_class = AuthorSerializer

# class BookList(generics.ListCreateAPIView):
#     queryset = Book.objects.all()
#     serializer_class = BookSerializer

# class CategoryList(generics.ListCreateAPIView):
#     queryset = Category.objects.all()
#     serializer_class = CategorySerializer
from django.shortcuts import get_object_or_404
from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status
from .models import Author, Book, Publisher

from django.shortcuts import get_object_or_404
from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status
from .models import Author, Book, Publisher

class AuthorAPIView(APIView):

    def get(self, request, pk=None):
        if pk is not None:
            author = get_object_or_404(Author, pk=pk)
            data = {
                'id': author.id,
                'name': author.name,
                'email': author.email,
                'books': [book.id for book in author.written_books.all()]
            }
            return Response(data, status=status.HTTP_200_OK)
        else:
            authors = Author.objects.all()
            data = []
            for author in authors:
                data.append({
                    'id': author.id,
                    'name': author.name,
                    'email': author.email,
                    'books': [book.id for book in author.written_books.all()]
                })
            return Response(data, status=status.HTTP_200_OK)

    def post(self, request):
        name = request.data.get('name')
        email = request.data.get('email')
        books = request.data.get('books')

        author = Author.objects.create(name=name, email=email)
        if books is not None:
            author.author_books.set(books)

        data = {
            'id': author.id,
            'name': author.name,
            'email': author.email,
            'books': [book.id for book in author.written_books.all()]
        }

        return Response(data, status=status.HTTP_201_CREATED)

class BookAPIView(APIView):

    def get(self, request, pk=None):
        if pk is not None:
            book = get_object_or_404(Book, pk=pk)
            data = {
                'id': book.id,
                'title': book.title,
                'description': book.description,
                'publisher': book.publisher.id,
                #'authors': [author.id for author in book.written_books.all()]
            }
            return Response(data, status=status.HTTP_200_OK)
        else:
            books = Book.objects.all()
            bookss = Author.objects.get(id=2)
            x=bookss.book_set.all()
            for i in x:
                print(i.publisher,bookss.name, i.description)
            data = []
            

            
            for book in books:
                pu=book.publisher
                print(pu.name)
                data.append({
                    'id': book.id,
                    'title': book.title,
                    'description': book.description,
                    'publisher': book.publisher.name,
                    'authors': [author.name for author in book.authors.all()]
                })
            return Response(data, status=status.HTTP_200_OK)

    def post(self, request):
        title = request.data.get('title')
        description = request.data.get('description')
        publisher_id = request.data.get('publisher')
        author_ids = request.data.get('authors')

        publisher = get_object_or_404(Publisher, pk=publisher_id)
        authors = Author.objects.filter(pk__in=author_ids)

        book = Book.objects.create(title=title, description=description, publisher=publisher)
        book.book_authors.set(authors)

        data = {
            'id': book.id,
            'title': book.title,
            'description': book.description,
            'publisher': book.publisher.id,
            'authors': [author.id for author in book.book_authors.all()]
        }

        return Response(data, status=status.HTTP_201_CREATED)


class PublisherAPIView(APIView):

    def get(self, request, pk=None):
        if pk is not None:
            publisher = get_object_or_404(Publisher, pk=pk)
            data = {
                'id': publisher.id,
                'name': publisher.name,
                'website': publisher.website,
                #'books': [book.id for book in publisher.books.all()]
            }
            return Response(data, status=status.HTTP_200_OK)
        else:
            publishers = Publisher.objects.all()
            data = []
            for publisher in publishers:
                data.append({
                    'id': publisher.id,
                    'name': publisher.name,
                    'website': publisher.website,
                    #'books': [book.id for book in publisher.books.all()]
                })
            return Response(data, status=status.HTTP_200_OK)

    def post(self, request):
        prin
        name = request.data.get('name')
        website = request.data.get('website')

        publisher = Publisher.objects.create(name=name, website=website)

        data = {
            'id': publisher.id,
            'name': publisher.name,
            'website': publisher.website,
            #'books': [book.id for book in publisher.books.all()]
        }

        return Response(data, status=status.HTTP_201_CREATED)


from rest_framework.views import APIView
from rest_framework.response import Response
from .models import Artist, Album, Song

class MusicView(APIView):
    def get(self, request):
        # Get all songs with their related artist and album models
        songs = Song.objects.select_related('artist', 'album').all()
        # Create a list of dictionaries to represent the song data
        song_data = []
        for song in songs:
            song_dict = {
                'id': song.id,
                'title': song.title,
                'artist': {
                    'id': song.artist.id,
                    'name': song.artist.name,
                    # Add other fields as needed
                },
                'album': {
                    'id': song.album.id,
                    'title': song.album.title,
                    # Add other fields as needed
                },
                # Add other fields as needed
            }
            song_data.append(song_dict)
        # Return the data as JSON
        return Response(song_data)
