# Generated by Django 4.0.4 on 2023-02-17 09:46

from django.db import migrations, models


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='userDetail',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(blank=True, max_length=10, null=True)),
                ('email_id', models.EmailField(blank=True, max_length=254, null=True)),
                ('phone_no', models.CharField(blank=True, max_length=10, null=True)),
                ('address', models.CharField(blank=True, max_length=50, null=True)),
            ],
        ),
    ]