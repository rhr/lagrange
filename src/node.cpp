/*
 * node.cpp
 *
 *  Created on: Nov 24, 2009
 *      Author: smitty
 */

#include <string>
#include <vector>
#include <cstring>
#include <iostream>
#include <sstream>
#include <map>

using namespace std;

#include "node.h"
#include "BranchSegment.h"
#include "node_object.h"
#include "string_node_object.h"
#include "vector_node_object.h"
#include "superdouble.h"

Node::Node():BL(0.0),height(0.0),number(0), name(""),
		parent(NULL),children(vector<Node *> ()),assoc(map<string,NodeObject *>()),
		assocDV(map<string,vector<Superdouble> >()),comment(""){
		}

Node::Node(Node * inparent):BL(0.0),height(0.0),number(0), name(""),
		parent(inparent),children(vector<Node *> ()),assoc(map<string,NodeObject *>()),
		assocDV(map<string,vector<Superdouble> >()),comment(""){
		}

Node::Node(double bl,int innumber,string inname, Node * inparent) :BL(bl),height(0.0),
		number(innumber), name(inname),parent(inparent),children(vector<Node *> ()),
		assoc(map<string,NodeObject *>()),assocDV(map<string,vector<Superdouble> >()),comment(""){
		}

vector<Node*> Node::getChildren(){
	return children;
}

bool Node::isExternal(){
	if(children.size()<1)
		return true;
	else
		return false;

}

bool Node::isInternal(){
	if(children.size()>0)
		return true;
	else
		return false;
}

bool Node::isRoot(){
	if(parent == NULL)
		return true;
	else
		return false;
}

bool Node::hasParent(){
	if(parent == NULL)
		return false;
	else
		return true;
}

void Node::setParent(Node & p){
	parent = &p;
}

int Node::getNumber(){
	return number;
}

void Node::setNumber(int n){
	number = n;
}

double Node::getBL(){
	return BL;
}

void Node::setBL(double bl){
	BL = bl;
}

double Node::getHeight(){
	return height;
}

void Node::setHeight(double he){
	height = he;
}

bool Node::hasChild(Node & test){
	bool ret = false;
	for(unsigned int i=0;i<children.size();i++){
		if(children.at(i)==&test){
			ret = true;
			break;
		}
	}
	return ret;
}

bool Node::addChild(Node & c){
	if(hasChild(c)==false){
		children.push_back(&c);
		c.setParent(*this);
		return true;
	}else{
		return false;
	}
}

bool Node::removeChild(Node & c){
	if(hasChild(c)==true){
		for(unsigned int i=0;i<children.size();i++){
			if(children.at(i)==&c){
				children.erase(children.begin()+i);
				break;
			}
		}
		return true;
	}else{
		return false;
	}
}

Node & Node::getChild(int c){
	return *children.at(c);
}

string Node::getName(){
	return name;
}

void Node::setName(string s){
	name = s;
}

string Node::getComment(){
	return comment;
}

void Node::setComment(string s){
	comment = s;
}

string Node::getNewick(bool bl){
	string ret = "";
	for(int i=0;i<this->getChildCount();i++){
		if(i==0)
			ret = ret+"(";
		ret = ret+this->getChild(i).getNewick(bl);
		if(bl==true){
			std::ostringstream o;
			o << this->getChild(i).getBL();
			ret = ret +":"+o.str();
		}
		if(i == this->getChildCount()-1)
			ret =ret +")";
		else
			ret = ret+",";
	}
	if(name.size()>0)
		ret = ret + name;
	return ret;
}

/*
 * should be returning the stringnodeobjects as the names for internal
 * nodes
 */
string Node::getNewick(bool bl,string obj){
	string ret = "";
	for(int i=0;i<this->getChildCount();i++){
		if(i==0)
			ret = ret+"(";
		ret = ret+this->getChild(i).getNewick(bl,obj);
		if(bl==true){
			std::ostringstream o;
			o << this->getChild(i).getBL();
			ret = ret +":"+o.str();
		}
		if(i == this->getChildCount()-1)
			ret =ret +")";
		else
			ret = ret+",";
	}
	if(isInternal()==true){
		if (obj == "number"){
			std::ostringstream o;
			o << number;
			ret = ret +o.str();
		}else{
			if(this->getObject(obj)!=NULL){
				std::ostringstream o;
				o << (*((StringNodeObject*) (this->getObject(obj))));
				ret = ret +o.str();
			}
		}
	}else{//EXTERNAL
		if(name.size()>0)
			ret = ret + name;
	}
	return ret;
}

/*
 * should return with branch lengths determined from the string obj
 */
string Node::getNewickOBL(string obj){
	string ret = "";
	for(int i=0;i<this->getChildCount();i++){
		if(i==0)
			ret = ret+"(";
		ret = ret+this->getChild(i).getNewickOBL(obj);
		std::ostringstream o;
		double bl = (*(VectorNodeObject<double>*) this->getChild(i).getObject(obj))[0];
		o << bl;
		ret = ret +":"+o.str();

		if(i == this->getChildCount()-1)
			ret =ret +")";
		else
			ret = ret+",";
	}
	if(name.size()>0)
		ret = ret + name;
	return ret;
}


Node * Node::getParent(){
	return parent;
}

int Node::getChildCount(){
	return children.size();
}

void Node::assocObject(string name,NodeObject & obj){
	//need to see why this doesn't work
  //if(assoc.count(name) > 0){ //take this out if doesn't work
  //delete assoc[name];
  //}
	assoc[name] = obj.clone();
}

void Node::assocDoubleVector(string name, vector<Superdouble> & obj){
	if (assocDV.count(name) > 0 ){
		assocDV.erase(name);
	}
	vector<Superdouble> tvec (obj.size());
	for (unsigned int i=0;i<obj.size();i++){
		tvec[i] = obj[i];
	}
	assocDV[name] = tvec;
}

vector<Superdouble> * Node::getDoubleVector(string name){
	return &assocDV[name];
}

void Node::deleteDoubleVector(string name){
	if (assocDV.count(name) > 0 ){
		assocDV.erase(name);
	}
}


void Node::initSegVector(){
	segs = new vector<BranchSegment> ();
}

vector<BranchSegment> * Node::getSegVector(){
	return segs;
}

void Node::deleteSegVector(){
	delete segs;
}

void Node::initExclDistVector(){
	excluded_dists = new vector<vector<int > > ();
}

vector<vector< int> > * Node::getExclDistVector(){
	return excluded_dists;
}

void Node::deleteExclDistVector(){
	delete excluded_dists;
}

/*
 * use the string ones like this
 * StringNodeObject sno("...a node object");
 * tree.getRoot()->assocObject("test",sno);
 * cout << *((StringNodeObject*) (tree.getRoot()->getObject("test"))) << endl;
 *
 * and the vector like
 * VectorNodeObject<int> vno;
 * vno.push_back(1);vno.push_back(2);
 * tree.getRoot()->assocObject("testvno",vno);
 * cout << ((VectorNodeObject<int> *) (tree.getRoot()->getObject("testvno")))->at(0) << endl;
 */

NodeObject  * Node::getObject(string name){
	return assoc[name];
}

