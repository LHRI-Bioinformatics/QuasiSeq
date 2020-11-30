from django.conf.urls import url, include
from django.contrib import admin
from django.conf.urls import patterns
from django.conf import settings
from quasiseq import views

urlpatterns = patterns( "" ,
        url(r'^$', views.index, name='index'),
        url(r'^fluData/', views.fluData, name = "fluData"),
#        url(r'^hivData/', views.hivData, name = "hivData"),
#        url(r'^Pro4mix/', views.Pro4mix, name = "Pro4mix"),
        url(r'^result/(?P<taskId>[-\w]+)/(?P<sampleName>[-\w]+)$', views.result, name='result'),
        url(r'^downloadFile/(?P<taskId>[-\w]+)/(?P<sampleName>[-\w]+)$', views.downloadFile, name='downloadFile'),

        #url(r'process/', views.process, name= "process"),

)
