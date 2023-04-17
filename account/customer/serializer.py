# from rest_framework import serializers
# from .models import Author, Book, Category

# class AuthorSerializer(serializers.ModelSerializer):
#     class Meta:
#         model = Author
#         fields = '__all__'

# class BookSerializer(serializers.ModelSerializer):
#     author = AuthorSerializer()
#     class Meta:
#         model = Book
#         fields = '__all__'

# class CategorySerializer(serializers.ModelSerializer):
#     books = BookSerializer(many=True)
#     class Meta:
#         model = Category
#         fields = '__all__'